// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#ifndef _VECGEN_OPT_HPP
#define _VECGEN_OPT_HPP

#include <map>
#include <set>
#include <list>
#include "vecgen.hpp"

class vecgen_opt {
public:
    vecgen_opt(vecgen *vg) : code(vg->code), vg(vg) {}
    void optimize();

private:
    typedef std::list<insn>::iterator iptr;
    typedef std::map<int,std::set<int> >::iterator phiptr;

    void to_ssa();
    void from_ssa();
    void dead_code();
    void copy_propagate();
    void loop_hoist();
    void cse();
    void immediate_fold();
    void mux_chains();

    bool is_output(int r);
    void insert_before(iptr ip, insn i) { code.insert(ip, i); mod = true; }
    void insert_after(iptr ip, insn i) { code.insert(++ip, i); mod = true; }
    void remove(iptr ip) { insn i = {}; *ip = i; mod = true; }
    void replace(int keep, int toss, bool phis_too=true);
    void replace_all(int keep, int toss);
    bool modcheck(int r, std::map<int,bool>& modified);

    int gen_phi(iptr i, int r, std::map<int,int>& dupmap);
    void get_regs_used(int r, std::map<int,bool>& used);

    std::map<int, std::set<int> > phi;    // phi id -> set of regs
    std::map<std::set<int>, int> phi_ids; // back mapping
    std::list<insn>& code;
    std::map<int,int> output_phis;      // output register to phi id
    std::map<int,int> output_reg_alias; // scratch reg to output reg
    vecgen* vg;
    bool mod;
};

#endif // _VECGEN_OPT_HPP
