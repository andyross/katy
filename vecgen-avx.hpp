// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#ifndef _VECGEN_AVX_HPP
#define _VECGEN_AVX_HPP

#include <map>
#include <list>

#include "avxgen.hpp"

// AVX-specific code generation.  Used by vecgen after the
// architecture-independent stages.  Uses an avxgen to generate the
// actual instructions.
//
// FIXME: the register allocator is pretty much entirely
// platform-neutral beyond an NREGS value.  Move somewhere else?
class vecgen_avx {
public:
    vecgen_avx(vecgen *vg) : vg(vg) {}
    void codegen();

private:
    enum mop { OP, SPILL, FILL, POOLTEST };

    // AVX code generation keeps an augmented instruction record with
    // physical registers, spill/fill marks, and pointers between loop
    // endpoints.
    struct avxinsn {
        avxinsn(mop o) : mo(o), pdst(-1), pa(-1), pb(-1), pc(-1) { insn i={}; vi=i; }
        avxinsn(insn i) : vi(i), mo(OP), pdst(-1), pa(-1), pb(-1), pc(-1) {}
        insn vi;
        mop mo;
        int pdst, pa, pb, pc;

        std::map<int,int> sphill;  // spill/fill virt->phys registers

        int addr; // generated code address
        std::list<avxinsn>::iterator loop_match;

        // Shorthand for SPILL/FILL, which reference just one reg
        int sphill_v() { return sphill.begin()->first; }
        int sphill_p() { return sphill.begin()->second; }
    };

    typedef std::list<avxinsn>::iterator aiptr;

    struct regs { int r[16]; };

    void dump();
    void dump_insn(avxinsn& i);
    void copy_code();
    void fixup_bools();
    void register_assignment();
    void loop_reset(aiptr i, regs curr, regs &old);
    void assign_scratch_slots();
    avxgen::RegMem locate_reg(avxinsn &i, int r, bool &simd);
    bool is_mutable(int r);
    int fill_reg(int vr, aiptr ip, regs& regs);
    int find_reg(int vr, aiptr ip, regs& regs);
    void prune_regs(aiptr ip, regs& regs, std::list<int> *rorder=0);
    void fill(aiptr ip, int vr, int r);
    void spill(aiptr ip, int vr, int r);
    void emit_load_store(avxgen& v, avxinsn& i);
    void emit_cull(avxgen& v, aiptr i);
    void emit_call(avxgen &v, int mi_cb, int mi_data);
    std::vector<int> cull_regs(struct insn &i);
    bool unassigned(aiptr ip, int r);
    bool unused(aiptr ip, int r);

    class vecgen *vg;
    std::list<avxinsn> code;
    std::map<int,int> scratch; // vr -> scratch slot
    int bool_true; // immediate register id for 0x80000000
    int nscratch; // max number of scratch registers on the stack
};

#endif // _VECGEN_AVX_HPP
