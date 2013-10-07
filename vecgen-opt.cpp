// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include <cstdio> // DEBUG
#include <stdexcept>
#include <stack>
#include <cmath>
#include <cstring>
#include "vecgen-opt.hpp"

using namespace std;

typedef list<insn>::iterator iptr;

static iptr loop_end(iptr i)
{
    if(i->o == LOOP) i++; // FIXME: dumb edge case due to last-- in gen_phi()
    for(int depth=0; depth > 0 || i->o != POOL; i++) {
        if(i->o == LOOP) depth++;
        if(i->o == POOL) depth--;
    }
    return i;
}

bool vecgen_opt::is_output(int r)
{
    return r == vg->msk || vg->outputs.hasid(r) || contains(output_reg_alias, r);
}

// Finds the sets of all possible assignments to "r" (detected via
// dupmap[dst]==r) from which this instruction can be reached.  Walk
// backwards to the top of the innermost enclosing loop: any
// assignment to r here MUST have been executed and thus any earlier
// ones would have been overriden, so return.  Reaching the top of the
// loop, start over at the end and walk back: the last assignment to r
// here MAY have been executed so add it to the set and continue.
// Then start over with the next enclosing loop until the begining of
// code has been reached.
int vecgen_opt::gen_phi(iptr i, int r, map<int,int>& dupmap)
{
    iptr last=i;
    last--; // FIXME/tricky: *i must be inspected because it can read itself
    int depth = 0;
    bool at_top = true;
    set<int> regs;
    while(i-- != code.begin()) {
        if(i->o == POOL) {
            depth++;
        } else if(i->o == LOOP && depth > 0) {
            depth--;
        } else if(i->o == LOOP && depth == 0) {
            iptr loopend = loop_end(last);
            for(iptr j = loopend; j != last; j--) {
                if(dupmap[j->dst] == r) {
                    regs.insert(j->dst);
                    break;
                }
            }
            last = ++loopend;
        } else if(dupmap[i->dst] == r) {
            regs.insert(i->dst);
            at_top = false;
            break;
        }
    }

    // If we reach the top of code looking, add the original
    // (input/output/constant/immediate) register.
    if(at_top)
        regs.insert(r);

    if(regs.size() == 1)
        return *regs.begin();

    if(contains(phi_ids, regs))
        return phi_ids[regs];

    int p = vg->nextid();
    phi[p] = regs;
    phi_ids[regs] = p;

    return p;
}

void vecgen_opt::to_ssa()
{
    // Renumber all assignments to be unique, tracking duplicate
    // assignments to the same initial register in dupmap[].  Outputs
    // already have one assignment before code starts.
    map<int,bool> seen;
    for(map<int,int>::iterator i=vg->outputs.i2v.begin(); i!=vg->outputs.i2v.end(); i++)
        seen[i->first] = true;
    seen[vg->msk] = true;

    map<int,int> dupmap; // dst register to "original" register
    for(iptr i=code.begin(); i != code.end(); i++) {
        if(i->dst && seen[i->dst]) {
            int id = vg->nextid();
            dupmap[id] = dupmap[i->dst] = i->dst;
            i->dst = id;
        }
        seen[i->dst] = true;
    }

    // Compute phi sets for all uses of duplicate registers
    for(iptr i=code.begin(); i != code.end(); i++) {
        if(i->a && dupmap[i->a])
            i->a = gen_phi(i, dupmap[i->a], dupmap);
        if(i->b && dupmap[i->b])
            i->b = gen_phi(i, dupmap[i->b], dupmap);
        if(i->c && dupmap[i->c])
            i->c = gen_phi(i, dupmap[i->c], dupmap);
    }

    // The output registers are "read" after the end of code.  We
    // aren't actually generating phis: the structure of "no forward
    // branch" code is such that only one final assignment to a
    // register can exist; gen_phi() is just an easy way of getting
    // that assignment.
    for(map<int,int>::iterator i=vg->outputs.i2v.begin(); i!=vg->outputs.i2v.end(); i++)
        if(dupmap[i->first])
            output_reg_alias[gen_phi(code.end(), dupmap[i->first], dupmap)] = i->first;
    if(dupmap[vg->msk])
        output_reg_alias[gen_phi(code.end(), dupmap[vg->msk], dupmap)] = vg->msk;
}

void vecgen_opt::from_ssa()
{
    // Build a reverse phi table that maps a register to a set of phi
    // ids that contain the register.
    std::map<int, std::set<int> > r2phi;
    for(phiptr i=phi.begin(); i!=phi.end(); i++)
        for(set<int>::iterator j=i->second.begin(); j!=i->second.end(); j++)
            r2phi[*j].insert(i->first);

    // For each assignment, check for any phis that reference it and
    // emit MOVs to do the copy.
    for(iptr i=code.begin(); i!=code.end(); i++) {
        if(i->dst && contains(r2phi, i->dst)) {
            set<int>& s = r2phi.find(i->dst)->second;
            for(set<int>::iterator j = s.begin(); j != s.end(); j++) {
                insn mov = { MOV, *j, i->dst };
                insert_after(i, mov);
            }
        }
    }

    // Emit copies to the beginning of code as phi "assignments".
    // OPTIMIZATION FIXME: this emits needless MOVs, should renumber
    // one of the assignments instead.
    for(int i=0; i<vg->next_id; i++) {
        if(contains(r2phi, i) && !vg->is_scratch(i)) {
            set<int>& s = r2phi.find(i)->second;
            for(set<int>::iterator j = s.begin(); j != s.end(); j++) {
                insn mov = { MOV, *j, i };
                insert_before(code.begin(), mov);
            }
        }
    }

    // Renumber the known output register aliases to the output registers
    map<int,int> ora2 = output_reg_alias;
    for(map<int,int>::iterator i=ora2.begin(); i!=ora2.end(); i++)
        replace_all(i->second, i->first);
}

void vecgen_opt::get_regs_used(int r, map<int,bool>& used)
{
    used[r] = true;
    if(contains(phi, r))
        for(set<int>::iterator i=phi[r].begin(); i!=phi[r].end(); i++)
            used[*i] = true;
}

void vecgen_opt::dead_code()
{
    map<int,bool> used, assigned;
    for(int i=0; i<vg->next_id; i++)
        used[i] = is_output(i);
    for(iptr i=code.begin(); i!=code.end(); i++) {
        if(vg->is_control(i->o))
            continue;
        assigned[i->dst] = true;
        get_regs_used(i->a, used);
        get_regs_used(i->b, used);
        get_regs_used(i->c, used);
    }
    for(map<int,int>::iterator i=output_phis.begin(); i!=output_phis.end(); i++)
        used[i->first] = true;

    // Prune unused instructions
    for(iptr i=code.begin(); i!=code.end(); i++)
        if(i->dst && !used[i->dst])
            remove(i);

    // Prune unused phis ids too
    list<int> deadphis;
    for(phiptr i=phi.begin(); i!=phi.end(); i++)
        if(!used[i->first])
            deadphis.push_back(i->first);
    for(list<int>::iterator i=deadphis.begin(); i!=deadphis.end(); i++)
        phi.erase(*i);
}

// Replaces references in operand fields to "toss" with "keep"
// "phis_too" indicates that we should also replace references to
// "toss" in (non-degenerate) phi sets (typical, but copy propagation
// doesn't want to do this).
void vecgen_opt::replace(int keep, int toss, bool phis_too)
{
    for(iptr i=code.begin(); i!=code.end(); i++) {
        if(i->a == toss) i->a = keep;
        if(i->b == toss) i->b = keep;
        if(i->c == toss) i->c = keep;
    }
    for(phiptr i=phi.begin(); i!=phi.end(); i++) {
        int id = i->first;
        if(contains(phi[id], toss)) {
            if(phi[id].size() > 1 && !phis_too)
                continue;
            set<int> regs = phi[id];
            regs.erase(toss);
            regs.insert(keep);
            if(contains(phi_ids, regs)) {
                // Already exists, do a recursive replace.
                replace(phi_ids[regs], id);
            } else {
                phi_ids.erase(regs);
                phi[id] = regs;
                phi_ids[regs] = id;
                mod = true;
            }
        }
    }
    if(contains(output_reg_alias, toss)) {
        output_reg_alias[keep] = output_reg_alias[toss];
        output_reg_alias.erase(toss);
    }
    if(contains(phi, toss))
        phi[toss].clear();
}

// Replaces dst references too.
void vecgen_opt::replace_all(int keep, int toss)
{
    replace(keep, toss);
    for(iptr i=code.begin(); i!=code.end(); i++) {
        if(i->dst == toss) {
            i->dst = keep;
            mod = true;
        }
    }
}

void vecgen_opt::copy_propagate()
{
    for(iptr i=code.begin(); i!=code.end(); i++)
        if(i->o == MOV && !is_output(i->dst))
            replace(i->a, i->dst, false);
}

// Is r (or its phi expansion) in modified[]?
bool vecgen_opt::modcheck(int r, map<int,bool>& modified)
{
    if(!contains(phi, r))
        return modified[r];
    for(set<int>::iterator i=phi[r].begin(); i!=phi[r].end(); i++)
        if(modified[*i])
            return true;
    return false;
}

// Walk forward making a stack of LOOP starts, and at each POOL pop
// one off the stack and test for assignments to values unread within
// the loop.  This gets us inside-to-outside-then-first-to-last order.
void vecgen_opt::loop_hoist()
{
    // Need a set of registers that appear in phi sets.  We can't
    // hoist phi assignments because the transitive copies end up in
    // the wrong place.  FIXME: we could still hoist the generator
    // instruction and add a copy here, or even put in a placeholder
    // "PHICOPY" instruction (which would also be useful for copy
    // propagation?).  Probably not worth it except for hoisting
    // thigns like RCP or RSQ...
    map<int,bool> phiregs;
    for(phiptr i=phi.begin(); i!=phi.end(); i++)
        for(set<int>::iterator j=i->second.begin(); j!=i->second.end(); j++)
            phiregs[*j] = true;

    stack<iptr> loops;
    for(iptr i=code.begin(); i!=code.end(); i++) {
        if(i->o == LOOP) {
            loops.push(i);
        } else if(i->o == POOL) {
            // Make a list of registers modified in the loop
            map<int,bool> mod;
            for(iptr j=loops.top(); j != i; j++)
                if(j->dst)
                    mod[j->dst] = true;

            // Test each instruction in the loop for dependence on the modified set
            for(iptr j=loops.top(); j != i; j++) {
                if(!j->a || j->o == POOL || j->o == IF || phiregs[j->dst])
                    continue;
                if(j->o == CULL || j->o == CULL_FLD)
                    continue; // culls can never hoist
                if(!modcheck(j->a, mod) && !modcheck(j->b, mod) && !modcheck(j->c, mod)) {
                    insert_before(loops.top(), *j);
                    remove(j);
                }
            }

            // If the POOL msk argument is *itself* invariant, that
            // means we have a loop that only executes once (or an
            // infinite loop, but we aren't responsible for correctly
            // generating non-terminating code!).  Remove the
            // LOOP/POOL entirely.  This also has the nice effect of
            // removing empty loops for free.
            if(!modcheck(i->a, mod)) {
                remove(loops.top());
                remove(i);
            }
            loops.pop();
        }
    }
}

// So insn can be a map key in cse()
static inline bool operator< (const insn& a, const insn& b)
{
    return memcmp(&a, &b, sizeof(a)) < 0;
}

void vecgen_opt::cse()
{
    map<insn,iptr> exprs;
    for(iptr i=code.begin(); i!=code.end(); i++) {
        if(i->dst && i->o != LOAD) {
            insn e = *i;
            e.dst = 0; // index by insn with zeroed dst
            if(contains(exprs, e))
                replace(exprs[e]->dst, i->dst);
            exprs[e] = i;
        }
    }
}

void vecgen_opt::immediate_fold()
{
    for(iptr i=code.begin(); i!=code.end(); i++) {
        int narg=0, nimm=0;
        narg += !!i->a;  nimm += !!vg->is_imm(i->a);
        narg += !!i->b;  nimm += !!vg->is_imm(i->b);
        narg += !!i->c;  nimm += !!vg->is_imm(i->c);
        if(!nimm)
            continue;

        if(nimm == narg) {
            // All operands are immediates
            if(i->o == MOV) continue; // noop optimization

            vecgen::uif a, b, c, dst;
            a.i = i->a ? vg->imm_val(i->a) : 0;
            b.i = i->b ? vg->imm_val(i->b) : 0;
            c.i = i->c ? vg->imm_val(i->c) : 0;
            dst.i = 0;

            // If it's interpretable, swap the instruction for a simple
            // MOV of a new immediate which will (probably) be further
            // optimized via copy propagation.
            if(vg->interpret_op(i->o, dst, a, b, c)) {
                insn mov = { MOV, i->dst, vg->rid(vg->imm_i(dst.i)) };
                *i = mov;
                mod = true;
            }
        }
    }
}

void vecgen_opt::mux_chains()
{
    map<int,iptr> muxes; // MUX operations by dst
    map<int, vector<iptr> > muxparents; // muxes that are parents/users of the key
    for(iptr i=code.begin(); i!=code.end(); i++) {
        if(i->o == MUX) {
            muxes[i->dst] = i;
            if(i->a) muxparents[i->a].push_back(i);
            if(i->b) muxparents[i->b].push_back(i);
            if(i->c) muxparents[i->c].push_back(i);
        }
    }

    // Chained MUX operators with the same selector operand are
    // degenerate: That is: a?b:(a?c:d) => a?b:d, and a?(a?:b:c):d =>
    // a?b:d
    for(map<int,iptr>::iterator mi=muxes.begin(); mi!=muxes.end(); mi++) {
        iptr i = mi->second;
        if(contains(muxes, i->a) && muxes[i->a]->c == i->c)
            i->a = muxes[i->a]->a;
        if(contains(muxes, i->b) && muxes[i->b]->c == i->c)
            i->b = muxes[i->b]->b;
    }

    // Similarly, if a "parent" (i.e. user) of an operation a MUX, and
    // an operand "child" is a MUX with the same selector, then the
    // child is degenerate and the parent can use a duplicated
    // operation without it (but leave the original in place as it
    // might be used elsewhere).  FIXME: not phi-aware, potentially
    // adds new instructions when the original cannot be removed.
    // FIXME: the mutation here can "break" a mux chain leaving
    // needless muxes that can't be removed.  Build a list of
    // modifications and apply them all at once?  Still needs
    // attention as to order...
    iptr i = code.end(); i--; do {
        for(unsigned int j=0; j<muxparents[i->dst].size(); j++) {
            iptr par = muxparents[i->dst][j];
            int sel = par->c;
            insn dup = *i;
            int *args[] = { &dup.a, &dup.b, &dup.c };
            for(int k=0; par->dst != i->dst && k<3; k++) {
                int *arg = args[k];
                if(*arg == par->dst || *arg == i->dst) {
                    continue; // cycle: don't touch
                } if(contains(muxes, *arg) && muxes[*arg]->c == sel) {
                    if(par->a == i->dst) {
                        par->a = dup.dst = vg->nextid();
                        *arg = muxes[*arg]->a;
                    } else if(par->b == i->dst) {
                        par->b = dup.dst = vg->nextid();
                        *arg = muxes[*arg]->b;
                    } else {
                        continue; // don't double-hit removed insns
                    }
                    insert_before(i, dup);
                    mod = true;
                }
            }
        }
    } while(--i != code.end());
}

void vecgen_opt::optimize()
{
    to_ssa();
    do {
        mod = false;
        immediate_fold();
        copy_propagate();
        loop_hoist();
        cse();
        mux_chains();
        dead_code();
    } while(mod);
    from_ssa();

    // Misc: yank noop A=MOV(A) instructions
    for(iptr i=code.begin(); i!=code.end(); i++)
        if(i->o == MOV && i->dst == i->a)
            remove(i);

    // FIXME: some elements of the immediates array may now be unused.
    // Need a pruning step.
}
