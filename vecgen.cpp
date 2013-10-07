// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include <algorithm>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "vecgen.hpp"
#include "vecgen-avx.hpp"
#include "vecgen-opt.hpp"

using namespace std;

typedef list<insn>::iterator iptr;

vr::vr(const vr& v)
{
    // G++ STL will copy-construct from default-constructed objects
    // (which sounds wrong to me) when initializing vectors, so check
    // for that as a synonym for "clone".
    if(v.cloneme || (v.mode == TMP && !v.rid)) {
        memcpy(this, &v, sizeof(v));
        cloneme = false;
    } else {
        init(REG);
        copyvg(&v);
        rid = vg->nextid();
        *this = v;
    }
}

vr& vr::operator=(const vr &ca)
{
    vr &a = const_cast<vr&>(ca); // we mutate temporaries
    if(a.mode == TMP && !a.rid)
        throw domain_error("assign from unused temporary");
    if(!vg && !a.vg)
        throw domain_error("vr assignment without relevant vecgen");
    copyvg(&a);
    if(mode == TMP && !rid)
        rid = vg->nextid();
    if(a.mode == OP)
        a.rid = rid;
    vg->gen(&a);
    if(a.mode != OP)
        vg->emit(MOV, rid, a.rid);
    return *this;
}

vr load(memsz sz, int obj, const vr &i)
{
    if(!i.vg)
        throw domain_error("load from anonymous index");
    return vrdup(vr(LOAD, sz, obj, &i.vg->mskreg, &i));
}

void store(memsz sz, int obj, const vr &i, const vr &v)
{
    vecgen *vg = i.vg ? i.vg : v.vg;
    if(!vg)
        throw domain_error("store() without vecgen pointer");
    vg->gen(&const_cast<vr&>(i));
    vg->gen(&const_cast<vr&>(v));
    vg->store(sz, obj, i.rid, v.rid);
}

vecgen::vecgen(const char *opts)
{
    if(opts)
        set_opts(opts);
    next_id = 1;
    num_culls = 0;
    msk = nextid();
    mskreg = vr(this, msk);
    flow.push_back(TOP);
    countidx = -1;
    nmems = 0;
}

void vecgen::gen(vr *v)
{
    if(!v || v->mode == vr::REG || (v->mode == vr::TMP && v->rid)) {
        return;
    } else if(v->mode == vr::IMM) {
        uif u;
        u.f = v->val;
        v->rid = gen_immediate(u.i);
        v->mode = vr::REG;
    } else if(v->mode == vr::IMMI) {
        v->rid = gen_immediate(v->ival);
        v->mode = vr::REG;
    } else if(v->mode == vr::OP) {
        if(!v->rid)
            v->rid = nextid();
        if(is_imm(v->rid)) throw domain_error("assign to immediate");
        if(inputs.hasid(v->rid)) throw domain_error("assign to input");
        if(consts.hasid(v->rid)) throw domain_error("assign to constant");
        gen(const_cast<vr*>(v->a));
        gen(const_cast<vr*>(v->b));
        gen(const_cast<vr*>(v->c));
        insn i = { v->o, v->rid, (v->a ? v->a->rid : 0), (v->b ? v->b->rid : 0),
                   (v->c ? v->c->rid : 0), v->imm0, v->imm1, 0 };
        emit(i);
    } else {
        throw domain_error("cannot generate temporary vr");
    }
}

// Trivial syntax: comma-separated tokens of the form "key=value" or
// just "key" (where value is inferred as "1")
void vecgen::set_opts(const char* optspec)
{
    string s = optspec;
    size_t start=0, end;
    do {
        end = s.find(',', start);
        string opt = s.substr(start, end == string::npos ? end : end-start);
        size_t eq = opt.find('=');
        string key = opt.substr(0, eq);
        string val = eq == string::npos ? "1" : opt.substr(eq+1);
        opts[key] = val;
        start = end+1;
    } while(end != string::npos);
}

const char* vecgen::opt(const char* key)
{
    return contains(opts, key) ? opts[key].c_str() : 0;
}

int vecgen::opt_i(const char* key)
{
    const char *v = opt(key);
    return v ? atoi(v) : 0;
}

bool vecgen::is_scratch(int r)
{
    return r != msk && !consts.hasid(r) && !inputs.hasid(r) && !outputs.hasid(r) && !is_imm(r);
}

void vecgen::store(memsz sz, int obj, int idx, int val)
{
    insn i = { STORE, 0, msk, idx, val, (int)sz, obj };
    emit(i);
}

int vecgen::cull(vr* regs, int nregs, int memidx_ptr, int memidx_count)
{
    // Implementation here needs some explanation: culls in the
    // instruction stream are one instruction per register, so that
    // the source argument tracking works in the optimizer.  The code
    // generator is expected to coalesce all sequential culls with the
    // same id and emit them along with the index increments at the
    // position of the first one (which gets a unique opcode).  One
    // optimizer gotcha: per-field culls can only hoist if *all* of
    // them hoist.
    int id = num_culls++;
    for(int i=0; i<nregs; i++) {
        insn o = { i ? CULL_FLD : CULL, 0, regs[i].rid, msk, 0, memidx_ptr, memidx_count, id };
        emit(o);
    }
    return id;
}

vr vecgen::imm(float f)
{
    uif u;
    u.f = f;
    return imm_i(u.i);
}

void vecgen::break_loop()
{
    if(find(flow.begin(), flow.end(), IN_LOOP) == flow.end())
        throw domain_error("break outside of loop");
    emit(BREAK);
}

void vecgen::end_loop()
{
    if(flow.back() != IN_LOOP)
        throw domain_error("end_loop outside of loop");
    // This turns into "jump-back-to-start-of-loop" at native codegen
    // time.  So it implicitly reads MSK to detect the termination
    // condition.  Make that source register explicit so that
    // optimization doesn't need special case handling.
    emit(POOL, 0, msk);
    flow.pop_back();
}

void vecgen::start_else()
{
    if(code.back().o == FI) {
        // ELSE right after FI: yank the FI and pretend we're still in the if
        flow.push_back(IN_IF);
        code.pop_back();
    } else if(flow.back() != IN_IF) {
        throw domain_error("start_else outside of if");
    }
    emit(ELSE);
    flow.pop_back();
    flow.push_back(IN_ELSE);
}

void vecgen::end_if()
{
    if(flow.back() != IN_IF && flow.back() != IN_ELSE)
        throw domain_error("end_if outside of if");
    if(flow.back() == IN_IF)
        emit(ELSE);
    emit(FI);
    flow.pop_back();
}

// STL has no predicate for "next/prev" that doesn't have side effects...
static iptr next(iptr i) { return ++i; }
static iptr prev(iptr i) { return --i; }

// Returns an iterator pointing to the end token of the nearest enclosing block
// (e.g. if you pass it an IF, it will return the matching FI).
static iptr block_end(iptr i)
{
    for(int depth=0; depth || (i->o != FI && i->o != POOL); i++) {
        if(i->o == IF || i->o == LOOP) depth++;
        if(i->o == FI || i->o == POOL) depth--;
    }
    return i;
}

// Note: no overrun checking, so use only if you know "o" exists in the list
static iptr find_at_level(iptr i, op o)
{
    for(int depth=0; depth || i->o != o; i++) {
        if(i->o == IF || i->o == LOOP) depth++;
        if(i->o == FI || i->o == POOL) depth--;
    }
    return i;
}

void vecgen::insert_before(iptr p, op o, int dst, int a, int b, int c)
{
    insn i = { o, dst, a, b, c };
    code.insert(p, i);
}

// Returns an iterator to the inserted element
iptr vecgen::insert_after(iptr p, op o, int dst, int a, int b, int c)
{
    insn i = { o, dst, a, b, c };
    code.insert(++p, i);
    return prev(p);
}

void vecgen::emit(insn &i)
{
    if(code.size() && code.back().o == BREAK)
        if(i.o != ELSE && i.o != FI && i.o != POOL)
            throw domain_error("break must be followed by end-of-block");
    code.insert(code.end(), i);
}

void vecgen::emit(op o, int dst, int a, int b, int c)
{
    insn i = { o, dst, a, b, c };
    emit(i);
}

// Brute force: predicate all instructions with a destination operand.
// The giant MUX nest will be unwound by the optimizer.
void vecgen::add_predication()
{
    map<int,bool> inuse;
    for(iptr i=code.begin(); i!=code.end(); i++) {
        if(i->dst && (inuse[i->dst] || is_output(i->dst))) {
            int d0 = i->dst, tmpdst = nextid();
            i->dst = tmpdst;
            i = insert_after(i, MUX, d0, d0, tmpdst, msk); // DST := MSK ? TMPDST : DST
        }
        inuse[i->dst] = true;
    }
}

// SIMD flow control is done via "emulated forward jumps".  There is a global
// predication/mask register "MSK" (register ID 1), which controls assignments
// whose results are needed outside of the block.  On reaching a point where a
// forward jump would normally be taken, the current value of MSK is saved off
// to a scratch register MSKn, MSK is set to 0, and on reaching the jump target
// MSKn is copied back to MSK (obviously all those actions are themselves
// predicated on MSK, allowing nested loops).  The one complexity is a BREAK,
// which is capable of "skipping over" the MSK restoration instructions of any
// IF/ELSE blocks between it and the enclosing loop.  So breaks also write
// zeros to the "saved msk" register of each skipped restore.  (Note,
// interestingly, that breaks don't need a save register themselves: if MSK is
// set, then they set it to zero.  If it's already zero, then obviously they
// won't affect it and it stays zero.  So the effect of BREAK is to clear MSK
// always).
void vecgen::add_forward_jumps()
{
    set<int> msksaves;
    for(iptr i = code.begin(); i != code.end(); i++) {
        if(i->o == IF) {
            // save & modify MSK at start, then insert copy back to MSK before the FI
            int predreg = i->a;
            int save = nextid();
            i = insert_after(i, MOV, save, msk);
            i = insert_after(i, AND, msk, save, predreg); // MSK = A & MSK

            // The ELSE starts over again from the saved mask and inverts the
            // modifier.
            iptr els = find_at_level(next(i), ELSE);
            if(next(els)->o != FI)
                insert_after(els, ANDN, msk, save, predreg); // MSK = SAVE & ~A

            // And restore msk
            insert_before(block_end(i), MOV, msk, save);
            msksaves.insert(save);
        } else if(i->o == LOOP) {
            // Loops are basically the same but simpler (there's no test at the
            // start, and breaks can't cross the end)
            int save = nextid();
            insert_before(i, MOV, save, msk);
            insert_after(block_end(next(i)), MOV, msk, save);
        } else if(i->o == BREAK) {
            // Find the save registers restorations that need to be cleared
            // (save = save ANDN msk, i.e. zero if msk is 1, else unchanged)
            for(iptr j=i; j->o != POOL; j++)
                if(j->o == MOV && j->dst == msk && contains(msksaves, j->a))
                    insert_before(i, ANDN, j->a, j->a, msk);
            insert_before(i, MOV, msk, imm_i(0).rid); // then clear msk
        }
    }
}

int vecgen::gen_immediate(int val)
{
    for(unsigned int i=0; i<immrecs.size(); i++)
        if(immrecs[i].val == val)
            return immrecs[i].rid;
    immrec ir = { nextid(), val };
    immrecs.push_back(ir);
    immr2i[ir.rid] = immrecs.size() - 1;
    return ir.rid;
}

string vecgen::regname(int r)
{
    ostringstream ss(ostringstream::out);
    ss <<   "R" << r;
    if     (r == msk)         ss << "/MSK";
    else if(is_imm(r))        ss << "/IMM" << immr2i[r];
    else if(outputs.hasid(r)) ss << "/OUT" << outputs.value(r);
    else if( consts.hasid(r)) ss <<   "/C" << consts.value(r);
    else if( inputs.hasid(r)) ss <<  "/IN" << inputs.value(r);
    return ss.str();
}

const char* vecgen::opname(op o)
{
    static map<op, const char*> opnames;
    static bool initialized = false;
    if(!initialized) {
#define O(o) opnames[o] = #o
        O(NOP); O(FADD); O(FSUB); O(FMUL); O(RCP); O(RSQ); O(MUX); O(FEQ);
        O(FNE); O(FLT); O(FGT); O(FLE); O(FGE); O(FLOOR); O(CEIL); O(F2I);
        O(I2F); O(AND); O(OR); O(ANDN); O(BITANDN); O(BITAND); O(BITOR);
        O(BITXOR); O(FMA); O(FMAX); O(FMIN); O(MOV); O(LOAD); O(STORE); O(IF);
        O(ELSE); O(FI); O(LOOP); O(BREAK); O(POOL); O(CULL); O(CULL_FLD);
        O(DEBUG_LOG);
#undef O
        initialized = true;
    }
    if(!opnames[o]) throw domain_error("unrecognized opcode");
    return opnames[o];
}

void vecgen::set_cull_bound(int cullid, int mi_max, int mi_data, int mi_callback)
{
    cullmeta cm = { mi_max, mi_data, mi_callback };
    cullbounds[cullid] = cm;
}

void vecgen::log(const char *msg, vr a, vr b)
{
    insn i = { DEBUG_LOG, 0, a.rid, b.rid, msk, (int)logmsgs.size() };
    emit(i);
    logmsgs.push_back(msg);
}

void vecgen::codegen()
{
    if(flow.size() != 1 || flow.back() != TOP)
        throw domain_error("codegen with unclosed blocks");

    add_predication();
    add_forward_jumps();

    if(!opt_i("no-optimize")) {
        vecgen_opt vo(this);
        vo.optimize();
    }

    vecgen_avx va(this);
    native_code.clear();
    va.codegen();
}

void vecgen::dump_code(FILE* log)
{
    fprintf(log, "[CODE]\n");
    for(iptr i = code.begin(); i != code.end(); i++) {
        if(i->o == NOP) continue;
        if(i->dst)
            fprintf(log, "  %-5s = ", regname(i->dst).c_str());
        else if(i->o == STORE || i->o == BREAK)
            fprintf(log, "  ");
        fprintf(log, "%s(", opname(i->o));
        if(i->a) fprintf(log, "%s", regname(i->a).c_str());
        if(i->b) fprintf(log, ", %s", regname(i->b).c_str());
        if(i->c) fprintf(log, ", %s", regname(i->c).c_str());
        fprintf(log, ")");
        if(i->imm0 || i->imm1 || i->imm2)
            fprintf(log, " [%d, %d, %d]", i->imm0, i->imm1, i->imm2);
        fprintf(log, "\n");
    }
    fprintf(log, "[END CODE]\n\n");
}

bool vecgen::interpret_op(op o, uif& dst, uif a, uif b, uif c)
{
    switch(o) {
    case FADD:     dst.f = a.f + b.f;                 break;
    case FSUB:     dst.f = a.f - b.f;                 break;
    case FMUL:     dst.f = a.f * b.f;                 break;
    case RCP:      dst.f = 1/a.f;                     break;
    case RSQ:      dst.f = 1/sqrt(a.f);               break;
    case AND:      dst.f = a.f != 0 && b.f != 0;      break;
    case OR:       dst.f = a.f != 0 || b.f != 0;      break;
    case MUX:      dst   = c.f == 0 ? a : b;          break;
    case FEQ:      dst.f = a.f == b.f;                break;
    case FNE:      dst.f = a.f != b.f;                break;
    case FLT:      dst.f = a.f < b.f;                 break;
    case FGT:      dst.f = a.f > b.f;                 break;
    case FLE:      dst.f = a.f <= b.f;                break;
    case FGE:      dst.f = a.f >= b.f;                break;
    case FLOOR:    dst.f = floor(a.f);                break;
    case CEIL:     dst.f = ceil(a.f);                 break;
    case F2I:      dst.i = (int)a.f;                  break;
    case I2F:      dst.f = (float)a.i;                break;
    case FMAX:     dst.f = a.f > b.f ? a.f : b.f;     break;
    case FMIN:     dst.f = a.f < b.f ? a.f : b.f;     break;
    case ANDN:     dst.f = (a.f != 0) && !(b.f != 0); break;
    case BITANDN:  dst.i = a.i & ~b.i;                break;
    case BITAND:   dst.i = a.i & b.i;                 break;
    case BITOR:    dst.i = a.i | b.i;                 break;
    case BITXOR:   dst.i = a.i ^ b.i;                 break;
    case FMA:      dst.f = a.f + b.f * c.f;           break;
    case MOV:      dst = a;                           break;
    case LOAD: case STORE: case IF: case ELSE: case FI:
    case LOOP: case POOL: case BREAK: case NOP:
    case CULL: case CULL_FLD: case DEBUG_LOG:
        return false;
    }
    return true;
}

void vecgen::interpret_cull(uif &msk, iptr ip, map<int,uif> &regs, void **mems, FILE *log)
{
    const int simd_width = 8; // FIXME: find a place to put this
    int *ptr = (int*)mems[ip->imm0];
    long idx = (long)mems[ip->imm1];
    if(msk.f != 0) {
        vector<int> cullregs;
        for(iptr j = code.begin(); j != code.end(); j++)
            if((j->o == CULL || j->o == CULL_FLD) && j->imm2 == ip->imm2)
                cullregs.push_back(j->a);
        for(unsigned int j=0; j<cullregs.size(); j++) {
            ptr[simd_width*j + idx] = regs[cullregs[j]].i;
            if(log) fprintf(log, "    cull field %d: %s=0x%08x{%g}\n",
                            j, regname(cullregs[j]).c_str(),
                            regs[cullregs[j]].i, regs[cullregs[j]].f);
        }
        if((idx = (idx+1) % simd_width) == 0)
            ptr += cullregs.size() * simd_width;
        mems[ip->imm0] = ptr;
        mems[ip->imm1] = (void*)idx;
        if(log) {
            fprintf(log, "    mem%d = %p (ptr)\n", ip->imm0, mems[ip->imm0]);
            fprintf(log, "    mem%d = %ld (count)\n", ip->imm1, (long)mems[ip->imm1]);
        }
    }
    if(contains(cullbounds, ip->imm2)) {
        cullmeta &c = cullbounds[ip->imm2];
        if(ptr == mems[c.mi_max] && idx) {
            typedef void(*cb_t)(void*);
            cb_t cb = reinterpret_cast<cb_t>(mems[c.mi_cb]);
            if(log) fprintf(log, "    Invoking cull bound callback.\n");
            cb(mems[c.mi_data]);
        }
    }
}

static const char* memsz_name(memsz m)
{
    switch(m) {
    case U8: return "U8";
    case S8: return "S8";
    case U16: return "U16";
    case S16: return "S16";
    case U32: return "U32";
    default: return "<unk>";
    }
}

void vecgen::interpret(float *constants, float *ins, float *outs, void **mems, FILE *log)
{
    map<int,uif> regs;
    regs[msk].f = 1;

    // Initialize the registers by ID
    for(vector<immrec>::iterator i=immrecs.begin(); i!=immrecs.end(); i++)
        regs[i->rid].i = i->val;
    for(map<int,int>::iterator i=consts.i2v.begin(); i!=consts.i2v.end(); i++)
        regs[i->first].f = constants[i->second];
    for(map<int,int>::iterator i=inputs.i2v.begin(); i!=inputs.i2v.end(); i++)
        regs[i->first].f = ins[i->second];

    stack<iptr> loops;
    for(iptr ip = code.begin(); ip != code.end(); ip++) {
        op o = ip->o;
        uif &dst=regs[ip->dst], a=regs[ip->a], b=regs[ip->b], c=regs[ip->c];

        if(o == DEBUG_LOG) {
            if(log) fprintf(log, "DEBUG_LOG %s, %s\n",
                            regname(ip->a).c_str(), regname(ip->b).c_str());
            if(c.f != 0) {
                // Note that debug logs always go out stdout, not the
                // register-level log stream which is usually disabled.
                printf("%s", logmsgs[ip->imm0].c_str());
                if(ip->a) printf(" 0x%x{%g}", a.i, a.f);
                if(ip->b) printf(" 0x%x{%g}", b.i, b.f);
                printf("\n");
            }
        } else if(o == NOP) {
            continue; // these are noops left by the optimizer
        } else if(is_control(ip->o)) {
            // These are noops during interpretation.  Log them anyway as landmarks.
            if(log) fprintf(log, "(%s)\n", opname(o));
        } else if(o == LOOP) {
            if(log) fprintf(log, "(LOOP)\n");
            loops.push(ip);
        } else if(o == POOL) {
            if(log) fprintf(log, "(POOL) %s = 0x%x{%g}\n",regname(ip->a).c_str(), a.i, a.f);
            if(a.f != 0)
                ip = loops.top();
            else {
                if(log) fprintf(log, "     (exit loop)\n");
                loops.pop();
            }
        } else if(o == LOAD) {
            memsz sz = (memsz)ip->imm0;
            int obj = ip->imm1;
            if(a.f) { // msk is in A
                if(sz == U8)  dst.i = ((int8_t*)mems[obj])[b.i] & 0xffu;
                if(sz == S8)  dst.i = ((int8_t*)mems[obj])[b.i];
                if(sz == U16) dst.i = ((int16_t*)mems[obj])[b.i] & 0xffffu;
                if(sz == S16) dst.i = ((int16_t*)mems[obj])[b.i];
                if(sz == U32) dst.i = ((int32_t*)mems[obj])[b.i];
                if(log) fprintf(log, "LOAD:%s %-5s = mem%d[%s=%d]\n           = %08x{%g}\n",
                                memsz_name(sz), regname(ip->dst).c_str(), obj,
                                regname(ip->b).c_str(), b.i, dst.i, dst.f);
            } else {
                if(log) fprintf(log, "LOAD:%s %-5s = mem%d[%s=%d] (masked: %s)\n",
                                memsz_name(sz), regname(ip->dst).c_str(), obj,
                                regname(ip->b).c_str(), b.i, regname(ip->a).c_str());
            }
        } else if(o == STORE) {
            memsz sz = (memsz)ip->imm0;
            int obj = ip->imm1;
            if(a.f) { // msk is in A
                if(sz == U8 || sz == S8)   ((int8_t*)mems[obj])[b.i] = c.i;
                if(sz == S16 || sz == S16) ((int16_t*)mems[obj])[b.i] = c.i;
                if(sz == U32)              ((int32_t*)mems[obj])[b.i] = c.i;
                if(log) fprintf(log, "STORE:%s mem%d[%d] = 0x%x{%g}\n",
                                memsz_name(sz), obj, b.i, c.i, c.f);
            } else {
                if(log) fprintf(log, "STORE:%s mem%d[%d] = 0x%x{%g} (masked: %s)\n",
                                memsz_name(sz), obj, b.i, c.i, c.f, regname(ip->a).c_str());
            }
        } else if(o == CULL) {
            if(log) fprintf(log, "CULL%d mem[%d,%d] msk=%g\n", ip->imm2, ip->imm0, ip->imm1, b.f);
            interpret_cull(b, ip, regs, mems, log);
        } else {
            // Traditional DST=OP(A,B,C) logic:
            if(is_imm(ip->dst) || inputs.hasid(ip->dst) || consts.hasid(ip->dst))
                throw domain_error("write to read-only register");

            interpret_op(o, dst, a, b, c);

            // Deliberately reduce precision for RCP/RSQ to approximate AVX behavior
            if(o == RCP || o == RSQ)
                dst.i &= ~0xfff;

            if(log) {
                fprintf(log, "%-5s = %s ( ", regname(ip->dst).c_str(), opname(o));
                if(ip->a) fprintf(log,   "%s=0x%08x{%g}", regname(ip->a).c_str(), a.i, a.f);
                if(ip->b) fprintf(log, ", %s=0x%08x{%g}", regname(ip->b).c_str(), b.i, b.f);
                if(ip->c) fprintf(log, ", %s=0x%08x{%g}", regname(ip->c).c_str(), c.i, c.f);
                fprintf(log, " )\n      = 0x%08x{%g}\n", dst.i, dst.f);
            }
        }
    }

    // Write the outputs back to the passed array
    for(map<int,int>::iterator i=outputs.i2v.begin(); i!=outputs.i2v.end(); i++)
        outs[i->second] = regs[i->first].f;
}
