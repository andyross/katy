// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include <map>
#include <set>
#include <list>
#include <stack>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>

#include "vecgen.hpp"
#include "avxgen.hpp"
#include "vecgen-avx.hpp"

using namespace std;

// Scalar registers used:
//
// RDI/R7  - First function argument: constants (scalar)
// RSI/R6  - Second function argument: inputs (simd)
// RDX/R2  - Third function argument: outputs (simd)
// RCX/R1  - Fourth function argument: array of memory pointers
// R8      - Fifth function argument: pointer to MSK array (simd)
// RSP/R4* - Stack pointer (conventional usage)
// RBP/R5* - Base pointer (conventional usage)
// R10     - Pointer to code-local immediates array
// RAX/R0  - Scalar index register for LOAD/STORE
//         - Iteration count for CULL
// RBX/R3* - "Memory" pointer for LOAD/STORE
//         - Temporary copy register for CULL
// R11     - Indexes pointer for LOAD/STORE
//         - Ouptut pointer for CULL
// R12*    - "Data" poitner for LOAD/STORE
//         - Output index for CULL
// R13*    - Scalar data register for LOAD/STORE
// R14*    - Mask pointer register for LOAD/STORE
//
// * These registers are callee-saved on the stack at function
//   entry/exit.
//
// Stack setup is conventional.  RBP is saved at function entry, saved
// values are pushed, the stack is aligned and then copied back to RBP
// which is used to index the SIMD "scratch" registers on the stack.
//
// OPTIMIZATION FIXME: two of the five saved registers could be
// eliminated by using RIP addressing for immediates instead of
// keeping a pointer in R10; and eliminating RBP in favor of the
// RSP-based addressing (which has slightly larger code size) for
// scratch registers.

enum { RAX=0, RCX=1, RDX=2,  RBX=3,  RSP=4,  RBP=5,  RSI=6,  RDI=7,
       R8=8,  R9=9,  R10=10, R11=11, R12=12, R13=13, R14=14, R15=15 };

const int NREGS = 16;

bool vecgen_avx::is_mutable(int r)
{
    return !vg->consts.hasid(r) && !vg->is_imm(r) && !vg->inputs.hasid(r);
}

void vecgen_avx::spill(aiptr ip, int vr, int r)
{
    // OPTIMIZATION FIXME: detect the case where we are spilling a
    // loop-invariant register and hoist the spill before the loop.
    // Note that this will break the existing scratch allocator
    // though!
    if(is_mutable(vr)) {
        avxinsn ai(SPILL);
        ai.sphill[vr] = r;
        code.insert(ip, ai);
    }
}

void vecgen_avx::fill(aiptr ip, int vr, int r)
{
    avxinsn ai(FILL);
    ai.vi.dst = vr;
    ai.sphill[vr] = r;
    code.insert(ip, ai);
}

// Detect registers that are unread in the future and zero them out.
// Optionally append to a list of the "read order" of used registers.
void vecgen_avx::prune_regs(aiptr ip, regs& regs, list<int> *rorder)
{
#define MARK(r) if(r && !assigned[r] && !read[r]) \
                { if(rorder) rorder->push_back(r); read[r] = true; }
    map<int,bool> assigned; // vregs definitely assigned in the future
    map<int,bool> read;     // vregs read in the future
    int depth = 0;
    for(aiptr start=ip, end=ip; end != code.end(); end++) {
        if(end->mo != OP && end->mo != POOLTEST) {
            continue;
        } else if(end->vi.o == LOOP) {
            depth++;
        } else if(end->vi.o == POOL && depth != 0) {
            depth--;
        } else if(end->vi.o == POOL) {
            // Unexpected POOL?  Find matching LOOP, walk from there
            // to start, update start.
            aiptr newstart = end->loop_match;
            for(aiptr l = newstart; l != start; l++) {
                if(l->mo == OP) {
                    MARK(l->vi.a);
                    MARK(l->vi.b);
                    MARK(l->vi.c);
                }
            }
            start = newstart;
        } else {
            // Normal instruction (or pooltest) not before "ip".
            // These instructions are definitely executed.  If they
            // lie in the future, note the assignment (which
            // suppresses future reads from entering the "used
            // farthest in future" list).  If this is ip itself, mark
            // it as a "read".  That prevents the register form being
            // incorrectly spilled mid-instruction even if it is not
            // used for a long time.
            MARK(end->vi.a);
            MARK(end->vi.b);
            MARK(end->vi.c);
            if(end->vi.dst) {
                if(ip == end) { MARK(end->vi.dst); }
                else          { assigned[end->vi.dst] = true; }
            }
        }
    }

    // Outputs are always read
    MARK(vg->msk);
    for(map<int,int>::iterator i=vg->outputs.i2v.begin(); i!=vg->outputs.i2v.end(); i++)
        MARK(i->first);

    // Clear slots that are known to be free
    for(int i=0; i<NREGS; i++)
        if(!read[regs.r[i]])
            regs.r[i] = 0;
#undef MARK
}

// Returns true if r has definitely *not* been assigned since its last
// FILL (or the beginning of code).
bool vecgen_avx::unassigned(aiptr ip, int r)
{
    aiptr last = ip;
    int depth = 0;
    do {
        ip--;
        if(ip->mo == OP && ip->vi.o == LOOP) {
            if(depth == 0) {
                for(last++; last->vi.o != POOL; last++)
                    if(ip->mo == OP && ip->vi.dst == r)
                        return false; // "might" have been assigned
            } else {
                depth--;
            }
        } else if(ip->mo == OP && ip->vi.o == POOL) {
            depth++;
        } else if(ip->vi.dst == r) {
            return ip->mo == FILL;
        }
    } while(ip != code.begin());
    return true;
}

// Returns true if the value of r is never read after ip (until either
// the next definite assignment or the end of code)
bool vecgen_avx::unused(aiptr ip, int r)
{
    if(vg->is_output(r) || r == vg->msk)
        return false;

    if(ip->vi.o == POOL)
        ip--; // hack, to allow use from loop_reset()

    aiptr loopcurr = ip;
    for(ip++; ip != code.end(); ip++) {
        if(ip->vi.a == r || ip->vi.b == r || ip->vi.c == r)
            return false;
        else if(ip->vi.dst == r)
            return ip->mo != FILL;

        if(ip->vi.o == POOL) {
            aiptr looptop = ip->loop_match; // loop_top(loopcurr);
            for(aiptr j=looptop; j != loopcurr; j++)
                if(j->vi.a == r || j->vi.b == r || j->vi.c == r)
                    return false;
                else if(j->vi.dst == r && j->mo == FILL)
                    return false;
                else if(j->vi.dst == r)
                    break;
            loopcurr = --looptop;
        }
    }
    return true;
}

// Really simple register allocator just grabs an empty register,
// discovers an unused one, or spills the one used farthest in the
// future (assumes loops are taken).  The walk order is a little
// interesting: when we reach the end of an enclosed loop, we bounce
// back to the top of it and walk to the earliest instruction yet read
// before moving on.  This therefore treats instructions in the
// current loop as "always before" ones after (i.e. it predicts all
// branches to be taken backward).
int vecgen_avx::find_reg(int vr, aiptr ip, regs& regs)
{
    list<int> rorder;       // earliest-to-latest order of vreg reads

    // First, see if it's already there
    for(int i=0; i<NREGS; i++)
        if(regs.r[i] == vr)
            return i;

    // Check for an empty slot, prune, then check again
    for(int i=0; i<2; i++) {
        for(int j=0; j<NREGS; j++) {
            if(!regs.r[j]) {
                regs.r[j] = vr;
                return j;
            }
        }
        if(i==0)
            prune_regs(ip, regs, &rorder);
    }

    // We need to choose a register to spill.  Walk backwards in
    // rorder to find the first one we see in regs[] (i.e. the
    // current register used farthest in the future).
    // OPTIMIZATION FIXME: rorder() doesn't take into account the
    // register reset at POOL, which will read registers "before" they
    // are used.
    for(list<int>::reverse_iterator i=rorder.rbegin(); i != rorder.rend(); i++) {
        for(int j=0; j<NREGS; j++) {
            if(regs.r[j] == *i) {
                if(!unassigned(ip, regs.r[j]) && !unused(ip, regs.r[j]))
                    spill(ip, regs.r[j], j);
                regs.r[j] = vr;
                return j;
            }
        }
    }
    throw logic_error("no spillable register");
}

// Ensures that a register contains vr
int vecgen_avx::fill_reg(int vr, aiptr ip, regs& regs)
{
    if(!vr)
        return -1;
    for(int i=0; i<NREGS; i++)
        if(regs.r[i] == vr)
            return i;
    int r = find_reg(vr, ip, regs);
    fill(ip, vr, r);
    return r;
}

// Returns storage location for a virtual register, and an boolean
// indicating whether the register is a calar or simd pointer.
//
// FIXME: use REG_RIP here for immediates, but that requires some way
// of tracking fixup locations on instructions yet to be emitted...
avxgen::RegMem vecgen_avx::locate_reg(avxinsn &i, int r, bool &simd)
{
    (void)i; // warning fix, i is passed for (disbled) scratch slot map
    int ptr, idx;
    simd = true;
    if     (vg->consts.hasid(r))  { ptr = 7;  idx = vg->consts.value(r); simd = false; }
    else if(vg->inputs.hasid(r))  { ptr = 6;  idx = vg->inputs.value(r);  }
    else if(vg->outputs.hasid(r)) { ptr = 2;  idx = vg->outputs.value(r); }
    else if(vg->is_imm(r))        { ptr = 10; idx = vg->immr2i[r]; simd = false; }
    else if(r == vg->msk)         { ptr = 8;  idx = 0; }
    else                          { ptr = 5;  idx = -1 - scratch[r]; }
    return avxgen::RegMem().base(ptr).disp((simd ? 32 : 4) * idx);
}

void vecgen_avx::emit_load_store(avxgen& v, avxinsn& i)
{
    // OPTIMIZATION FIXME: there are instruction bytes to be saved:
    // arranging for all the operands to be in low registers would
    // save 4 bytes of REX prefixes; and the 4-byte Jcc immediate
    // could be replaced with the 1-byte form with some hackery to
    // avxgen.  That would reduce the 136-byte sequence to 81, which
    // seems significant.  But remember the SNB uOp cache is in front
    // of the decoders and will act to hide extra insruction byte
    // overhead.

    int sz = i.vi.imm0, obj = i.vi.imm1;
    bool simd_index, simd_data, dummy;
    int datareg = i.vi.o == LOAD ? i.vi.dst : i.vi.c;
    v.mov64(RBX, v.mem(RCX).disp(8*obj));          // memory[obj] pointer
    v.lea(R11, locate_reg(i, i.vi.b, simd_index)); // index pointer
    v.lea(R12, locate_reg(i, datareg, simd_data)); // data pointer
    v.lea(R14, locate_reg(i, i.vi.a, dummy));      // mask pointer

    // If we're storing a scalar, read it into R13 before the loop
    // so we can skip it inside.
    if(!simd_data)
        v.mov32(R13, v.mem(R12).scale(2)); // do scalar data above loop
    for(int j=0; j<8; j++) {
        // Test the high bit of msk and skip if it's unset
        // FIXME: this is a hard-coded MSK register, but the actual
        // mask at runtime can be a scratch register...
        v.cmp32_imm8(v.mem(R14).disp(j*4), 0);
        int jmp_fixup = v.jns(0);

        v.mov32(RAX, v.mem(R11).disp(simd_index ? j*4 : 0)); // load index
        if(i.vi.o == LOAD) {
            // Load the indexed data from memory:
            if(sz == U8)  v.movzx8 (R13, v.mem(RBX).idx(RAX).scale(0));
            if(sz == S8)  v.movsx8 (R13, v.mem(RBX).idx(RAX).scale(0));
            if(sz == U16) v.movzx16(R13, v.mem(RBX).idx(RAX).scale(1));
            if(sz == S16) v.movsx16(R13, v.mem(RBX).idx(RAX).scale(1));
            if(sz == U32) v.mov32  (R13, v.mem(RBX).idx(RAX).scale(2));
            // Store it back to our dst register:
            v.mov32_store(v.mem(R12).disp(j*4), R13);
        } else {
            // Load the data to store:
            if(simd_data)
                v.mov32(R13, v.mem(R12).disp(j*4).scale(2));
            // Write it to dst memory + index:
            if(sz == U8 || sz == S8)    v.mov8_store (v.mem(RBX).idx(RAX).scale(0), R13);
            if(sz == U16 || sz == S16 ) v.mov16_store(v.mem(RBX).idx(RAX).scale(1), R13);
            if(sz == U32)               v.mov32_store(v.mem(RBX).idx(RAX).scale(2), R13);
        }
        v.fixup(jmp_fixup, v.size() - (jmp_fixup + 4));
    }
}

vector<int> vecgen_avx::cull_regs(struct insn &i)
{
    vector<int> cr;
    for(aiptr ip=code.begin(); ip!=code.end(); ip++)
        if(ip->mo == OP
           && (ip->vi.o == CULL || ip->vi.o == CULL_FLD)
           && ip->vi.imm2 == i.imm2)
            cr.push_back(ip->vi.a);
    return cr;
}

void vecgen_avx::emit_call(avxgen &v, int mi_cb, int mi_data)
{
    // All ymm registers are caller-save.  It would be better to spill
    // these to their proper storage instead of making more space on
    // the stack.  But we don't have access to register state here.
    v.mov64(RAX, RSP);      // save stack pointer
    v.and64_imm8(RSP, ~31); // and align
    v.sub64_imm32(RSP, 8*4*NREGS);
    for(int i=0; i<NREGS; i++)
        v.movaps_store(v.mem(RSP).disp(8*4*i), i);
    v.push(RAX); // put pre-alignment RSP on the stack below the regs

    // Caller-save registers: RCX, RDX, RSI, and R8. (RAX and R9-11
    // are also clobbered, but we don't use them in contexts where a
    // call can be made).  Argument goes in RDI.
    v.push(RCX);
    v.push(RDX);
    v.push(RSI);
    v.push(R8);
    v.push(RDI);
    v.mov64(RDI, v.mem(RCX).disp(8*mi_data));
    v.call(v.mem(RCX).disp(8*mi_cb));
    v.pop(RDI);
    v.pop(R8);
    v.pop(RSI);
    v.pop(RDX);
    v.pop(RCX);

    v.pop(RAX);
    for(int i=0; i<NREGS; i++)
        v.movaps(i, v.mem(RSP).disp(8*4*i));
    v.mov64(RSP, RAX);
}

// FIXME: need to be sure all input registers are spilled
void vecgen_avx::emit_cull(avxgen& v, aiptr i)
{
    vector<int> regs = cull_regs(i->vi);

    v.mov64(R11, v.mem(RCX).disp(8*i->vi.imm0)); // output pointer
    v.mov64(R12, v.mem(RCX).disp(8*i->vi.imm1)); // output count
    v.mov64_imm32(RAX, 0); // simd iteration counter
    int looptop = v.size();
    {
        bool simd;
        // Test msk and skip
        avxgen::RegMem rm = locate_reg(*i, i->vi.b, simd);
        v.cmp32_imm8(rm.idx(RAX).scale(2), 0);
        int jmp_fixup = v.jns(0);
        {
            // Copy all the fields
            for(unsigned int j=0; j<regs.size(); j++) {
                rm = locate_reg(*i, regs[j], simd);
                v.mov32(RBX, simd ? rm.idx(RAX).scale(2) : rm);
                v.mov32_store(v.mem(R11).disp(32*j).idx(R12).scale(2), RBX);
            }
            // Bump output counter, check if full, bump pointer if so
            v.inc(R12);
            v.cmp32_imm8(R12, 8);
            int jmp2_fixup = v.jnz(0);
            {
                v.mov64_imm32(R12, 0);
                v.add64_imm32(R11, 32*regs.size());
            }
            v.fixup(jmp2_fixup, v.size() - (jmp2_fixup+4));
        }
        v.fixup(jmp_fixup, v.size() - (jmp_fixup+4));
    }
    // Bump index, test for exit, repeat from beginning
    v.inc(RAX);
    v.cmp32_imm8(RAX, 8);
    v.jnz(looptop);

    // Save output pointer and count
    v.mov64_store(v.mem(RCX).disp(8*i->vi.imm0), R11);
    v.mov64_store(v.mem(RCX).disp(8*i->vi.imm1), R12);

    if(contains(vg->cullbounds, i->vi.imm2)) {
        vecgen::cullmeta &cm = vg->cullbounds[i->vi.imm2];

        // If we are at the maximum and have stored some data, then fire
        // the callback: we might overflow next time.
        v.cmp64(R11, v.mem(RCX).disp(8*cm.mi_max));
        int jmp1 = v.jnz(0);
        v.cmp32_imm8(R12, 0);
        int jmp2 = v.jz(0);
        emit_call(v, cm.mi_cb, cm.mi_data);
        v.fixup(jmp1, v.size() - (jmp1+4));
        v.fixup(jmp2, v.size() - (jmp2+4));
    }
}

// "Operand type" abstraction used by fixup_bools().  Note use of
// bitwise relationship: BOTH == BOOL | FLOAT
enum { NONE=0, BOOL=1, FLOAT=2, BOTH=3 };
static void optype(op o, int& dst, int& a, int& b, int& c)
{
    // FIXME: MOV and MUX "copy" their operands
    dst = a = b = c = FLOAT;
    switch(o) {
    case MUX:       c = BOOL; break;
    case POOL:      a = BOOL; break;
    case CULL:      b = BOOL; break;
    case LOAD:      a = BOOL; break;
    case STORE:     a = BOOL; break;
    case DEBUG_LOG: c = BOOL; break;
    case AND: case OR: case ANDN:
        dst = a = b = BOOL;
        break;
    case FNE: case FEQ: case FLT: case FLE: case FGT: case FGE:
        dst = BOOL;
        break;
    default:
        break;
    }
}

static bool ignore(op op)
{
    return op == IF || op == ELSE || op == FI || op == BREAK;
}

void vecgen_avx::fixup_bools()
{
    // OPTIMIZATION FIXME: float_zero can more efficiently be computed
    // as an XORPS of any register with itself.  Convert at FILL time,
    // and skip the allocation of the zero here?
    int float_zero = vg->rid(vg->imm(0));
    int float_one = vg->rid(vg->imm(1));

    // Build a table of registers read/written as bool/float and
    // decide on how that value will be "stored"
    map<int,int> read_type, write_type, stored_type;
    for(aiptr i=code.begin(); i!=code.end(); i++) {
        if(i->mo == OP) {
            if(ignore(i->vi.o))
                continue;

            int dt, at, bt, ct;
            optype(i->vi.o, dt, at, bt, ct);
            write_type[i->vi.dst] |= dt;
            read_type[i->vi.a] |= at;
            read_type[i->vi.b] |= bt;
            read_type[i->vi.c] |= ct;
        }
    }

    // Build a stored_type value for every register
    for(int r=1; r<vg->next_id; r++) {
        int rt = read_type[r], wt = write_type[r];
        if     (rt && rt != BOTH) stored_type[r] = rt;
        else if(wt && wt != BOTH) stored_type[r] = wt;
        else                      stored_type[r] = FLOAT;
    }

    // External values (except msk) are always FLOAT, the BOOL
    // representation is AVX-only
    for(map<int,int>::iterator i=stored_type.begin(); i!=stored_type.end(); i++)
        if(!vg->is_scratch(i->first))
            stored_type[i->first] = FLOAT;
    stored_type[vg->msk] = BOOL;

    for(aiptr i=code.begin(); i!=code.end(); i++) {
        if(i->mo != OP || ignore(i->vi.o))
            continue;

        int dt, at, bt, ct;
        optype(i->vi.o, dt, at, bt, ct);

        // Convert args on input (loop over a/b/c to avoid duplicated code)
        int *regs[] = { &i->vi.a, &i->vi.b, &i->vi.c };
        int types[] = { at, bt, ct };
        for(int j=0; j<3; j++) {
            int *r = regs[j], want = types[j];
            if(*r && want && want != stored_type[*r]) {
                insn fne = { FNE, 0, *r, float_zero };
                insn mux = { MUX, 0, float_zero, float_one, *r };
                insn &cvt = want == BOOL ? fne : mux;
                *r = cvt.dst = vg->nextid();
                code.insert(i, avxinsn(cvt));
            }
        }

        // Convert destination register after instruction
        aiptr next = i; next++;
        if(i->vi.dst && dt == BOOL && stored_type[i->vi.dst] == FLOAT) {
            // bool to float: DST = MUX(0, 1, DST)
            insn mux = { MUX, i->vi.dst, float_zero, float_one, i->vi.dst };
            code.insert(next, avxinsn(mux));
        } else if(i->vi.dst && dt == FLOAT && stored_type[i->vi.dst] == BOOL) {
            // float to bool: DST = FNE(DST, 0)
            insn fne = { FNE, i->vi.dst, i->vi.dst, float_zero };
            code.insert(next, avxinsn(fne));
        }

        // Skip over any instructions we might have added after i
        i = next; i--;
    }
}

static string prspec(int r)
{
    stringstream ss;
    if(r >= 0)
        ss << "@ymm" << r;
    return ss.str();
}

void vecgen_avx::dump_insn(avxinsn& i)
{
    //cout << hex << "{" << i.addr << "} " << dec;
    if(i.mo == SPILL || i.mo == FILL)
        cout << ((i.mo == FILL) ? "FILL " : "SPILL ")
             << vg->regname(i.sphill_v()) << " "
             << ((i.mo == FILL) ? "->"  : "<-") << " ymm" << i.sphill_p() << endl;
    else if(i.mo == POOLTEST)
        cout << "POOLTEST " << vg->regname(i.vi.a) << prspec(i.pa)
             << " " << vg->regname(i.vi.b) << prspec(i.pb) << endl;
    if(i.mo != OP)
        return;

    if(i.vi.dst) cout << vg->regname(i.vi.dst) << prspec(i.pdst) << " = ";
    cout << vecgen::opname(i.vi.o);
    if(i.vi.o == LOAD || i.vi.o == STORE || i.vi.imm0 || i.vi.imm1 || i.vi.imm2)
        cout << " [" << i.vi.imm0 << ", " << i.vi.imm1 << ", " << i.vi.imm2 << "] ";
    cout << "(";
    if(i.vi.o == IF || i.vi.o == POOL) {
        cout << vg->regname(i.vi.a);
    } else if(i.vi.o == CULL || i.vi.o == CULL_FLD) {
        cout << vg->regname(i.vi.a) << ", " << vg->regname(i.vi.b);
    } else {
        if(i.vi.a) cout << vg->regname(i.vi.a) << prspec(i.pa);
        if(i.vi.b) cout << ", " << vg->regname(i.vi.b) << prspec(i.pb);
        if(i.vi.c) cout << ", " << vg->regname(i.vi.c) << prspec(i.pc);
    }
    cout << ")";
    cout << endl;
}

void vecgen_avx::dump()
{
    // Assumes the presence of a objdump binary and a writable /tmp,
    // but should fail silently if those aren't there and just dump
    // the avxinsn code without interleaved disassembly.
    ofstream("/tmp/vecgen-avx.bin", ofstream::binary)
        .write((char*)&vg->native_code[0], vg->code_size());
    system("objdump -b binary -M x86-64 -m i386 -D "
           " /tmp/vecgen-avx.bin > /tmp/vecgen-avx.asm");

    aiptr ip = code.begin();
    ifstream assm("/tmp/vecgen-avx.asm");
    while(assm.good()) {
        size_t sep;
        int addr;
        string l;
        getline(assm, l);
        if((sep = l.find(":")) != string::npos) {
            istringstream(l.substr(0, sep)) >> hex >> addr;
            for(/**/; ip != code.end() && ip->addr <= addr; ip++)
                if(ip->mo != OP || ip->vi.o != NOP)
                    dump_insn(*ip);
        }
        cout << l << "\n";
    }
    for(/**/; ip != code.end(); ip++)
        dump_insn(*ip);
}

void vecgen_avx::copy_code()
{
    stack<aiptr> loops;
    for(vecgen::iptr i=vg->code.begin(); i!=vg->code.end(); i++) {
        if(ignore(i->o)) {
            code.push_back(avxinsn(*i)); // noop, but leave for debugging
        } else if(i->o == FMA) {
            // Convert FMA to FMUL/FADD. Two variants: one works in
            // place if dst and a are distict, the other needs a
            // scratch register.
            insn add = {}, mul = {};
            mul.o = FMUL;
            add.o = FADD;
            if(i->dst != i->a) {
                mul.dst = i->dst;  mul.a = i->b;    mul.b = i->c;
                add.dst = i->dst;  add.a = i->dst;  add.b = i->a;
            } else {
                mul.dst = vg->nextid(); mul.a = i->b;    mul.b = i->c;
                add.dst = i->dst;       add.a = mul.dst; add.b = i->a;
            }
            code.push_back(mul);
            code.push_back(add);
        } else if(i->o == LOOP) {
            code.push_back(avxinsn(*i));
            loops.push(--code.end());
        } else if(i->o == POOL) {
            // End-of-loop actually produces two instructions.  The
            // first represents the vcmpps that tests msk against
            // bool_true.  After that comes the register reset and
            // jnz back to loop top represented by the second
            // instruction.
            code.push_back(avxinsn(POOLTEST));
            code.back().vi.a = i->a;
            code.back().vi.b = bool_true;
            code.push_back(avxinsn(*i));
            aiptr pool = --code.end();
            loops.top()->loop_match = pool;
            pool->loop_match = loops.top();
            loops.pop();
        } else {
            code.push_back(avxinsn(*i));
        }
    }
}

// Resets register state at the end of the loop.  Updates the
// registers in "old" to reflect any changes.
void vecgen_avx::loop_reset(aiptr ip, regs curr, regs &old)
{
    // Emit spills for any registers whose contents need to be modified
    for(int i=0; i<NREGS; i++)
        if(old.r[i] && curr.r[i] && old.r[i] != curr.r[i])
            spill(ip, curr.r[i], i);


    // Copy from existing registers where possible
    bool clobbered[NREGS] = {};
    for(int i=0; i<NREGS; i++) {
        for(int j=0; j<NREGS; j++) {
            if(old.r[i] && curr.r[i] && i != j && old.r[i] == curr.r[j] && !clobbered[j]) {
                insn mov = { MOV, old.r[i], curr.r[j] };
                avxinsn ai(mov);
                ai.pdst = i;
                ai.pa = j;
                code.insert(ip, ai);
                clobbered[i] = true;
            }
        }
    }

    // Fill whatever is left
    for(int i=0; i<NREGS; i++)
        if(old.r[i] && curr.r[i] && old.r[i] != curr.r[i] && !clobbered[i])
            if(!unused(ip, old.r[i]))
                fill(ip, old.r[i], i);

    // Update "old" to reflect any new registers in curr
    for(int i=0; i<NREGS; i++)
        if(!old.r[i])
            old.r[i] = curr.r[i];
}

void vecgen_avx::register_assignment()
{
    // Need a table of cull fields by id
    map<int, vector<int> > culls;
    for(aiptr i=code.begin(); i!=code.end(); i++)
        if(i->mo == OP && (i->vi.o == CULL_FLD || i->vi.o == CULL))
            culls[i->vi.imm2].push_back(i->vi.a);

    stack<regs> regstack;
    regs zero_regs = {};
    regstack.push(zero_regs);
    for(aiptr i=code.begin(); i!=code.end(); i++) {
        regs &rtop = regstack.top();
        if(i->mo != OP && i->mo != POOLTEST) {
            continue;
        } else if(ignore(i->vi.o)) {
            // noop
        } else if(i->vi.o == LOOP) {
            // Dup the register state so we can restore at end of
            // loop.  Prune it first so we don't spend effort
            // restoring dead registers.
            prune_regs(i, regstack.top());
            regstack.push(regstack.top());
        } else if(i->vi.o == POOL) {
            // register state reset at end of loop
            regs curr = regstack.top();
            regstack.pop();
            loop_reset(i, curr, regstack.top());
        } else if(i->vi.o == CULL_FLD) {
            // noop
        } else if(i->vi.o == LOAD || i->vi.o == STORE) {
            // Inputs to these must be spilled if they're in
            // registers, outputs (from LOAD) appear in memory and
            // should be dropped if it's in a register.
            //
            // OPTIMIZATION FIXME: this loses track of the inputs and
            // forces them to be re-filled needlessly if they're used
            // again.  Find a way to note that the value is both in
            // memory and registers...
            for(int pr=0; pr<NREGS; pr++) {
                int vr = rtop.r[pr];
                if(vr && (vr == i->vi.a || vr == i->vi.b || vr == i->vi.c)) {
                    spill(i, vr, pr);
                    rtop.r[pr] = 0;
                }
                if(i->vi.o == LOAD && vr == i->vi.dst)
                    rtop.r[pr] = 0;
            }
        } else if(i->vi.o == CULL) {
            // Similar to STORE, all fields (including MSK in operand
            // B) must be purged from registers so they live in memory
            for(int pr=0; pr<NREGS; pr++) {
                int vr = rtop.r[pr];
                if(vr && vr == i->vi.b) {
                    spill(i, vr, pr);
                    rtop.r[pr] = 0;
                }
                for(unsigned int ci=0; ci<culls[i->vi.imm2].size(); ci++) {
                    if(vr && vr == culls[i->vi.imm2][ci]) {
                        spill(i, vr, pr);
                        rtop.r[pr] = 0;
                    }
                }
            }
        } else {
            // Regular instructions (and pooltest).
            i->pa = fill_reg(i->vi.a, i, rtop);
            i->pb = fill_reg(i->vi.b, i, rtop);
            i->pc = fill_reg(i->vi.c, i, rtop);
            if(i->vi.dst > 0)
                i->pdst = find_reg(i->vi.dst, i, rtop);
        }
    }

    // Add spills at the end of code to flush any output registers
    // left in regs[]
    regs &rt = regstack.top();
    for(int i=0; i<NREGS; i++)
        if(rt.r[i] == vg->msk || vg->outputs.hasid(rt.r[i]))
            spill(code.end(), rt.r[i], i);
}

// OPTIMIZATION FIXME: this was originally supposed to be an
// optimizing allocator that could reuse slots where the spills/fills
// don't overlap.  But it didn't work and was replaced with this 1:1
// hack.  Fix that, because for large functions there can be *many*
// very short-lived scratched registers and spreading them out on the
// stack hurts cache locality.
void vecgen_avx::assign_scratch_slots()
{
    nscratch = 0;
    for(aiptr i=code.begin(); i!=code.end(); i++) {
        if(i->mo == SPILL || i->mo == FILL) {
            int r = i->sphill_v();
            if(r && vg->is_scratch(r) && !contains(scratch, r))
                scratch[r] = nscratch++;
        }
    }
}

void vecgen_avx::codegen()
{
    bool_true = vg->rid(vg->imm_i(0x80000000));

    // Copy input code to an avxinsn list
    copy_code();

    // Booleans in the input representation are C-like, but in AVX
    // they are represented with the sign/high bit of the 32 bit word.
    // We need to add conversion logic when a boolean op takes a
    // non-boolean input, or vice versa.
    fixup_bools();

    // Do the register assignments, generating SPILL/FILL marks
    register_assignment();

    // Build a table mapping scratch registers to indexes so they can
    // be located on the stack.
    assign_scratch_slots();

    // OPTIMIZATION FIXME: add an optimization pass here where we
    // detect the case of a FILL of a register that is only used once.
    // If possible, fold that into the instruction as a direct memory
    // reference instead.

    // OPTIMIZATION FIXME: add a reordering pass here that sorts the
    // code by latency

    // Emit code!
    avxgen v(vg->native_code);
    if(vg->opt_i("break"))
        v.int3();
    v.push(RBP);     // push caller's RBP, then RBX, R12, R13
    v.push(RBX);
    v.push(R12);
    v.push(R13);
    v.push(R14);
    if(nscratch) {
        v.mov64(RAX, RSP);                // save off prealignment ("poppable") stack
        v.and64_imm8(RSP, ~31);           // align stack to 32 bytes
        v.mov64(RBP, RSP);                // establish stack frame
        v.sub64_imm32(RSP, 8*4*nscratch); // make room for scratch simd regs
        v.push(RAX);                      // push pre-alignment RSP
    } else {
        v.mov64(RBP, RSP);                // establish stack frame
    }

    v.lea(R10, v.mem().rip(0)); // immediates pointer into R10
    int imm_rip_fixup = v.size();

    int count_loop_top = v.size();
    stack<int> loop_tops;
    for(aiptr i=code.begin(); i!=code.end(); i++) {
        i->addr = vg->native_code.size();
        if((i->mo == SPILL || i->mo == FILL)) {
            bool simd;
            avxgen::RegMem rm = locate_reg(*i, i->sphill_v(), simd);
            if(i->mo == SPILL && is_mutable(i->sphill_v()))
                v.movaps_store(rm, i->sphill_p());
            else if(simd)
                v.movaps(i->sphill_p(), rm);
            else
                v.broadcastss(i->sphill_p(), rm);
        } else if(i->mo == POOLTEST) {
            v.testps(i->pa, i->pb); // msk vs. bool_true
        } else if(i->mo == OP) {
            switch(i->vi.o) {
            case LOOP:
                loop_tops.push(v.size());
                break;
            case POOL:
                // The flag setting happened earlier (see POOLTEST)
                v.jnz(loop_tops.top());
                loop_tops.pop();
                break;
            case LOAD: case STORE:
                emit_load_store(v, *i);
                break;
            case CULL:
                emit_cull(v, i);
                break;
            case FADD:  v.addps(i->pdst, i->pa, i->pb); break;
            case FSUB:  v.subps(i->pdst, i->pa, i->pb); break;
            case FMUL:  v.mulps(i->pdst, i->pa, i->pb); break;
            case RCP:   v.rcpps(i->pdst, i->pa); break;
            case RSQ:   v.rsqrtps(i->pdst, i->pa); break;
            case MUX:   v.blendvps(i->pdst, i->pa, i->pb, i->pc); break;
            case FEQ:   v.cmpps(i->pdst, i->pa, i->pb, avxgen::EQ); break;
            case FNE:   v.cmpps(i->pdst, i->pa, i->pb, avxgen::NEQ); break;
            case FLT:   v.cmpps(i->pdst, i->pa, i->pb, avxgen::LT); break;
            case FGT:   v.cmpps(i->pdst, i->pa, i->pb, avxgen::GT); break;
            case FLE:   v.cmpps(i->pdst, i->pa, i->pb, avxgen::LE); break;
            case FGE:   v.cmpps(i->pdst, i->pa, i->pb, avxgen::GE); break;
            case FLOOR: v.roundps(i->pdst, i->pa, avxgen::DOWN); break;
            case CEIL:  v.roundps(i->pdst, i->pa, avxgen::UP); break;
            case F2I:   v.cvttps2dq(i->pdst, i->pa); break;
            case I2F:   v.cvtdq2ps(i->pdst, i->pa); break;
            case FMAX:  v.maxps(i->pdst, i->pa, i->pb); break;
            case FMIN:  v.minps(i->pdst, i->pa, i->pb); break;
            case MOV:   v.movaps(i->pdst, i->pa); break;
            case BITXOR:   v.xorps(i->pdst, i->pa, i->pb); break;

            // OPTIMIZATION FIXME: float zero/non-zero booleans work
            // fine with bitwise AND/OR, so don't need coercing to AVX
            // booleans if it's just for these instructions.
            case AND:
            case BITAND:   v.andps(i->pdst, i->pa, i->pb); break;
            case OR:
            case BITOR:    v.orps(i->pdst, i->pa, i->pb); break;

            case ANDN:
            case BITANDN:  v.andnps(i->pdst, i->pb, i->pa); break;

            case FMA:
                throw domain_error("True FMA unimplemented");
            case NOP: case IF: case ELSE: case FI: case BREAK: case CULL_FLD:
            case DEBUG_LOG:
                break; // noop
            }
        }
    }

    // If we're in count mode, decrement and loop back.
    if(vg->countidx >= 0) {
        v.dec64(v.mem(RCX).disp(8*vg->countidx));
        int exit_fixup = v.jz(0);
        if(vg->instride)
            v.add64_imm32(RSI, 8*4*vg->instride);
        if(vg->outstride)
            v.add64_imm32(RDX, 8*4*vg->outstride);
        v.jmp(count_loop_top);
        v.fixup(exit_fixup, v.size() - (exit_fixup + 4));
    }

    if(nscratch) {
        v.pop(RAX);        // restore pre-alignment...
        v.mov64(RSP, RAX); // ... stack pointer
    }
    v.pop(R14);
    v.pop(R13);
    v.pop(R12);
    v.pop(RBX);
    v.pop(RBP); // pop caller's stack frame
    v.ret();

    // Now align to 4 bytes and postpend the scalar immediates
    while(v.size() & 3)
        v.data8(0x90);
    v.fixup(imm_rip_fixup - 4, v.size() - imm_rip_fixup);
    for(unsigned int i=0; i<vg->immrecs.size(); i++)
        v.data32(vg->immrecs[i].val);

    if(vg->opt_i("dump-asm"))
        dump();
}
