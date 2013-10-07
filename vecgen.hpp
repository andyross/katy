// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#ifndef _VECGEN_HPP
#define _VECGEN_HPP

//
// Top-level code generation for vector SIMD code.  A vecgen object
// can create vr ("virtual register") objects representing known
// storage and compose them together in C++ expressions that will
// automatically generate code to evaluate that expression.  The
// resulting code (after copying to executable memory) can be invoked
// through a C function pointer with known signature.
//

#include <list>
#include <string>
#include <vector>
#include <map>
#include <deque>
#include <cstring>

enum op { NOP=0, FADD, FSUB, FMUL, RCP, RSQ, MUX, FEQ, FNE, FLT, FGT, FLE, FGE,
          FLOOR, CEIL, F2I, I2F, AND, OR, ANDN, BITANDN, BITAND, BITOR, BITXOR, FMA, FMAX,
          FMIN, MOV, LOAD, STORE, IF, ELSE, FI, LOOP, BREAK, POOL, CULL, CULL_FLD, DEBUG_LOG };

enum memsz { U8, S8, U16, S16, U32 };

class vecgen;

// "Virtual register".  These are C++ objects that can be passed
// around and built into expressions using traditional syntax.  Those
// expressions automatically generate their code when used as an
// rvalue in an assignment or as an argument to start_if().
class vr {
public:
    // Default-constructed vr's become scratch registers on assignment
    vr() { init(TMP); }

    // Copy constructor produces a scratch register that gets assigned
    // from the rhs unless (!) the rhs is the output of one of vecgen
    // register factories (input/output/constant/immediate), which get
    // cloned.
    vr(const vr& v);

    // Converts from float to a runtime immediate.  Note: float, not
    // int. Use vecgen::imm_i to get a bit pattern.
    vr(float f) { init(IMM); val=f; }

    vr operator+(const vr &a)  const { return vrdup(vr(FADD,   this, &a)); }
    vr operator-(const vr &a)  const { return vrdup(vr(FSUB,   this, &a)); }
    vr operator*(const vr &a)  const { return vrdup(vr(FMUL,   this, &a)); }
    vr operator&(const vr &a)  const { return vrdup(vr(BITAND, this, &a)); }
    vr operator|(const vr &a)  const { return vrdup(vr(BITOR,  this, &a)); }
    vr operator&&(const vr &a) const { return vrdup(vr(AND,    this, &a)); }
    vr operator||(const vr &a) const { return vrdup(vr(OR,     this, &a)); }
    vr operator^(const vr &a)  const { return vrdup(vr(BITXOR, this, &a)); }
    vr operator==(const vr &a) const { return vrdup(vr(FEQ,    this, &a)); }
    vr operator!=(const vr &a) const { return vrdup(vr(FNE,    this, &a)); }
    vr operator<(const vr &a)  const { return vrdup(vr(FLT,    this, &a)); }
    vr operator<=(const vr &a) const { return vrdup(vr(FLE,    this, &a)); }
    vr operator>(const vr &a)  const { return vrdup(vr(FGT,    this, &a)); }
    vr operator>=(const vr &a) const { return vrdup(vr(FGE,    this, &a)); }

    // A precision warning about recip/rsqrt: on AVX, these can return
    // exact zeros for very large inputs (>= 2^126) instead of finite
    // values.
    friend vr recip(const vr &a)                  { return vrdup(vr(RCP,     &a)); }
    friend vr rsqrt(const vr &a)                  { return vrdup(vr(RSQ,     &a)); }
    friend vr andnot(const vr &a, const vr &b)    { return vrdup(vr(ANDN,    &a, &b)); }
    friend vr bitandnot(const vr &a, const vr &b) { return vrdup(vr(BITANDN, &a, &b)); }
    friend vr fma(const vr &a, const vr &b, vr c) { return vrdup(vr(FMA,     &a, &b, &c)); }
    friend vr floor(const vr &a)                  { return vrdup(vr(FLOOR,   &a)); }
    friend vr ceil(const vr &a)                   { return vrdup(vr(CEIL,    &a)); }
    friend vr f2i(const vr &a)                    { return vrdup(vr(F2I,     &a)); }
    friend vr i2f(const vr &a)                    { return vrdup(vr(I2F,     &a)); }
    friend vr max(const vr &a, const vr &b)       { return vrdup(vr(FMAX,    &a, &b)); }
    friend vr min(const vr &a, const vr &b)       { return vrdup(vr(FMIN,    &a, &b)); }
    friend vr mux(const vr &sel, const vr &a, const vr &b) { return vrdup(vr(MUX, &b, &a, &sel)); }

    vr& operator=(const vr &a);
    vr& operator+=(const vr &a) { return *this = *this + a; }
    vr& operator-=(const vr &a) { return *this = *this - a; }
    vr& operator*=(const vr &a) { return *this = *this * a; }
    vr& operator&=(const vr &a) { return *this = *this & a; }
    vr& operator|=(const vr &a) { return *this = *this | a; }
    vr& operator^=(const vr &a) { return *this = *this ^ a; }

    void setvg(vecgen *v) { vg = v; }

    // Memory operations.  Remember the i argument should be an
    // integer, e.g. from f2i().  Note also that because the integer
    // comes out of float math, it's limited to 24 bits of precision
    // (i.e. 16/32/64MB maximum for 1/2/4 byte elements) unless
    // special care is taken to composite integers via bitwise math.
    friend vr load(memsz sz, int obj, const vr &i);
    friend void store(memsz sz, int obj, const vr &i, const vr &val);

    friend vr imm(float f);
    friend vr imm_i(int i);

private:
    friend class vecgen;
    typedef enum { REG, OP, IMM, IMMI, TMP } vrmode;

    // This exists to defeat the return value optimization.  We want
    // "vr a = b+c;" to produce "a" as a scratch register (via the
    // copy constructor from "b+c"), not the FADD operation.  But RVO
    // is legally allowed to construct "b+c" directly in the memory
    // allocated for a, which is not what we want.
    friend vr vrdup(const vr &a) { return a; }

    vr(op op, const vr *ap=0, const vr *bp=0, const vr *cp=0)
    { init(OP); o=op; copyvg(ap, bp, cp); a=ap; b=bp; c=cp; }
    vr(op op, memsz sz, int obj, const vr *ap, const vr *bp=0)
    { init(OP); o=op; copyvg(ap, bp); imm0=sz; imm1=obj; a=ap; b=bp; }
    vr(vecgen* v, int r) { init(REG); vg=v; rid=r; cloneme=true; }

    void init(vrmode m) { memset(this, 0, sizeof(vr)); mode = m; }
    void copyvg(const vr *a=0, const vr *b=0, const vr *c=0)
    { if(!vg) { vg = a->vg ? a->vg : (b->vg ? b->vg : c->vg); } }

    vrmode mode;
    int rid;
    vecgen *vg;
    const vr *a, *b, *c;
    float val;
    int ival;
    op o;
    int imm0, imm1;
    bool cloneme;
};

inline vr imm(float f) { return vr(f); }
inline vr imm_i(int i) { vr v; v.mode = vr::IMMI; v.ival = i; return v; }

inline vr recip2(const vr &a) { vr b=recip(a); return (imm(2)-a*b)*b; }
inline vr rsqrt2(const vr &a) { vr b=rsqrt(a); return (b*0.5)*(imm(3)-a*b*b); }

template<class A, class B>
inline bool contains(const A &a, const B &b) { return a.find(b) != a.end(); }

// Trivial bidirectional id/value mapping for register files.
struct regmap {
    int probe(int value, int& nextid) {
        if(v2i.find(value) == v2i.end()) {
            v2i[value] = nextid++;
            i2v[v2i[value]] = value;
        }
        return v2i[value];
    }
    bool hasid(int id) { return i2v.find(id) != i2v.end(); }
    int value(int id) { return i2v[id]; }
    std::map<int,int> i2v, v2i;
};

struct insn {
    op o;
    int dst, a, b, c;
    int imm0, imm1, imm2;
};

typedef void (*vecgen_fn)(float* consts, float* ins, float* outs, void **mems, float *msk);

class vecgen {
public:
    vecgen(const char *opts=0);

    // Register generators.  Immediates are indexed and unified by
    // their values; inputs and outputs are indexed arrays passed to
    // the generated function.  Temporary scratch registers can be
    // created as needed, or will be created via the expression
    // generators.
    vr imm_i(int val) { return vr(this, gen_immediate(val)); }
    vr imm(float f);
    vr constant(int idx) { return vr(this, consts.probe(idx, next_id)); }
    vr input(int idx) { return vr(this, inputs.probe(idx, next_id)); }
    vr output(int idx) { return vr(this, outputs.probe(idx, next_id)); }
    vr scratch() { return vr(this, nextid()); }

    // Newton-Raphson refinements of the low-precision ops
    vr recip2(const vr &a) { vr b=recip(a); return (imm(2)-a*b)*b; }
    vr rsqrt2(const vr &a) { vr b=rsqrt(a); return (b*0.5)*(imm(3)-a*b*b); }

    void start_loop() { emit(LOOP); flow.push_back(IN_LOOP); }
    void break_loop();
    void end_loop();
    void start_if(vr &a) { gen(&a); emit(IF, 0, a.rid); flow.push_back(IN_IF); }
    void start_else();
    void end_if();

    // Simple sequential allocation of memory array (the "mems"
    // argument to the generated code) indexes.  Different code
    // generation modules need to be able to share this space without
    // clobbering each other.
    int alloc_memidx() { return nmems++; }
    int n_mems() { return nmems; }

    // Emit a scalar hsoa "slice" for each unmasked SIMD thread.
    // Takes indices into the mem[] array, one of which is an aligned
    // pointer to a presized array of hsoa records, the other is an
    // integer indicating how many entries [0:N-1] in the current
    // record are already filled.  Both are mutated by the generated
    // code and can thus be reused by another iteration.  The intent
    // of this API is to efficiently implement a "culling" step, where
    // one SIMD task (say, a vertex shader) conditionally generates
    // output records for another (a fragment shader).
    int cull(vr* regs, int nregs, int memidx_ptr, int memidx_count);

    // By default, cull regions are unbounded, they will continue
    // writing to memory as long as there is new data.  Bounds can be
    // set at any time before codegen().  This takes three indexes
    // into the mems[] array.  The first is a "maximum" pointer value
    // that the cull pointer is allowed to have.  After it is full the
    // function pointer stored at mi_callback (cdecl, signature
    // "void(*)(void*)") will be invoked with the uninspected pointer
    // value found at mi_cbdata.  The callback might be expected to
    // reset the ptr/count values with a new memory region, etc...
    void set_cull_bound(int cullid, int mi_max, int mi_cbdata, int mi_callback);

    // Insert a log message into the interpreter stream.  Ignored in
    // the generated code (except that it does affect optimizer
    // behavior as it reads registers).  Simply prints the message
    // followed by the values of the 0-2 arguments.
    void log(const char *msg, vr a=vr(), vr b=vr());

    // The returned code can be invoked through a standard cdecl
    // function pointer with the following signature:
    //
    //   void code(const float* consts, const simd_t* ins, simd_t* outs, void **mems, simd_t *msk)
    //
    // All "simd_t" pointers should be aligned (i.e. 32 bytes for
    // AVX).  The scalar and memory pointers need alignment only if
    // the platform requires it (x86 does not) or for performance
    // reasons.  Caller is responsible for sizing and allocating the
    // output buffer in executable memory.
    void codegen();
    int code_size() { return native_code.size(); }
    void get_code(char* buf) { std::memcpy(buf, &native_code[0], code_size()); }

    // The generated code can be configured to iterate over an array
    // of records.  The count is passed in at call time via the mems[]
    // array in the index specified.  The input and output base
    // pointers are incremented at each iteration by the specified
    // strides (in units of a SIMD register, e.g. 32 bytes).
    // FIXME: should this just be done via the opts string?
    void set_count(int count_idx, int input_stride, int output_stride)
    { countidx = count_idx; instride = input_stride; outstride = output_stride; }

    // Quick & dirty interpreter for the VM.  Makes debugging the
    // optimizer much easier and provides a good reference for test
    // suites.  Takes scalar (!) pointers: emulates logic, not SIMD.
    // Also doesn't support counted iteration.
    void interpret(float *consts, float *ins, float *outs, void **mems, FILE *log=0);

    void dump_code(FILE* log);

    // FIXME: used only by the unit test right now.  If we expose
    // these, it wants a better API than blindly expecting the user to
    // tightly pack.  Find the highest index and/or make the user set
    // it explicitly?
    int ninputs() { return inputs.i2v.size(); }
    int noutputs() { return outputs.i2v.size(); }

    bool is_immediate(const vr &v) { return is_imm(v.rid); }

private:
    friend class vr;
    friend class vecgen_avx;
    friend class vecgen_opt;
    friend vr load(memsz sz, int obj, const vr &i);
    friend void store(memsz sz, int obj, const vr &i, const vr &val);

    typedef std::list<insn>::iterator iptr;

    void emit(op o, int dst=0, int a=0, int b=0, int c=0);
    void emit(insn &i);

    int load(memsz sz, int obj, int idx);
    void store(memsz sz, int obj, int idx, int val);

    void set_opts(const char* opts);

    // Ensures the vr has a register ID and recursively emits code to
    // get it the right value at runtime.
    void gen(vr *v);
    int rid(const vr &v) { return v.rid; }
    int nextid() { return next_id++; }
    int gen_immediate(int val);

    void insert_before(iptr p, op o, int dst=0, int a=0, int b=0, int c=0);
    iptr insert_after(iptr p, op o, int dst=0, int a=0, int b=0, int c=0);

    void add_predication();
    void add_forward_jumps();

    union uif { int i; float f; };
    bool interpret_op(op o, uif& dst, uif a, uif b, uif c);
    void interpret_cull(uif &msk, iptr ip, std::map<int,uif> &regs, void **mems, FILE *log);

    bool is_scratch(int r);
    bool is_output(int r) { return outputs.hasid(r); }
    bool is_imm(int r) { return immr2i.find(r) != immr2i.end(); }
    int imm_val(int r) { return immrecs[immr2i[r]].val; }

    static bool is_control(op o) { return o==IF || o==ELSE || o==FI || o==BREAK; }

    const char* opt(const char* key);
    int opt_i(const char* key);
    std::string regname(int r);
    static const char* opname(op o);

    struct immrec { int rid, val; };
    enum flowstate { TOP, IN_IF, IN_ELSE, IN_LOOP };
    struct cullmeta { int mi_max, mi_data, mi_cb; };

    regmap inputs, outputs, consts;
    std::vector<immrec> immrecs;
    std::map<int,int> immr2i;
    std::deque<flowstate> flow;
    int next_id;
    int msk; // mask register
    vr mskreg;
    int num_culls;
    std::list<insn> code;
    std::vector<unsigned char> native_code;
    std::map<std::string,std::string> opts;
    std::map<int,cullmeta> cullbounds;
    int countidx;
    int instride, outstride;
    std::vector<std::string> logmsgs;
    int nmems;
};

// Helper macros that make block structure look more natural
inline int _startloop(vecgen& vg) { vg.start_loop(); return 1; }
inline int _endloop(vecgen& vg) { vg.end_loop(); return 0; }
inline int _startif(vecgen& vg, const vr &r) { vg.start_if(const_cast<vr&>(r)); return 1; }
inline int _startelse(vecgen& vg) { vg.start_else(); return 1; }
inline int _endif(vecgen& vg) { vg.end_if(); return 0; }
#define VRIF(vg, e) for(int _1st=_startif(vg, e); _1st; _1st=_endif(vg))
#define VRELSE(vg) for(int _1st=_startelse(vg); _1st; _1st=_endif(vg))
#define VRLOOP(vg) for(int _1st=_startloop(vg); _1st; _1st=_endloop(vg))

#endif // _VECGEN_HPP
