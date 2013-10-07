// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#ifndef _AVXGEN_HPP
#define _AVXGEN_HPP

#include <vector>
#include <string>

// avxgen: x64_64 code generator for VEX-encoded AVX/SSE instructions
// (plus a handful of scalar insructions vector code needs)

// Assembles x86-64 AVX instructions
const int REG_NONE = -1;
const int REG_RIP = -2;

class avxgen {
public:
    avxgen(std::vector<unsigned char> &c) : _code(c) {}

    void reset() { _code.clear(); }
    int size() { return _code.size(); }
    const unsigned char* code() { return &_code[0]; }

    // Test predicates.  There are more, but we don't bother with the
    // NaN result mapping or signalling.
    enum test { EQ=0, LT=1, LE=2, UNORD=3, NEQ=4, GE=5, GT=6, ORD=7 };

    // Rounding modes
    enum round { NEAR=0, DOWN=1, UP=2, TRUNC=3 };

    // Shorthand gadget for constructing x86 ModRM/SIB/disp
    // references.  Converts from an integer register specification as
    // you'd expect, or use the default constructor and chain methods
    // (base/idx/scale/disp) to set fields.  For example (given an
    // avxgen object named "v"), to increment the dst register by the
    // Nth element of a SIMD array pointed to by R8:
    //     v.addps(dst, dst, v.mem().base(8).idx(32*N));
    // Or to load a scalar constant from R9 into all SIMD slices of YMM2:
    //     v.broadcastss(2, v.mem().base(9).idx(constid).scale(2);
    // Or a RIP-relative(-to-the-next-insn) constant:
    //     v.broadcastss(3, v.mem().rip(const_offset));
    struct RegMem {
        char _reg, _base, _idx, _scale;
        int _disp;
        RegMem(int r) : _reg(r), _base(REG_NONE),_idx(REG_NONE), _scale(0), _disp(0) {}
        RegMem() : _reg(REG_NONE), _base(REG_NONE), _idx(REG_NONE), _scale(0), _disp(0) {}
        RegMem& base(char r)  { _base = r;  return *this; }
        RegMem& idx(char r)   { _idx = r;   return *this; }
        RegMem& scale(char s) { _scale = s; return *this; }
        RegMem& disp(int d)   { _disp = d;  return *this; }
        RegMem& rip(int d)  { _base = REG_RIP; _disp = d; return *this; }
    };

    // AVX instructions
    void addps(int dst, int src0, RegMem src1)
    { vex_nds(0x58, MM_0F, W_IGNORE, L_256, PP_NONE, dst, src0, src1); }
    void subps(int dst, int src0, RegMem src1)
    { vex_nds(0x5c, MM_0F, W_IGNORE, L_256, PP_NONE, dst, src0, src1); }
    void mulps(int dst, int src0, RegMem src1)
    { vex_nds(0x59, MM_0F, W_IGNORE, L_256, PP_NONE, dst, src0, src1); }
    void rcpps(int dst, RegMem src)
    { vex_nds(0x53, MM_0F, W_IGNORE, L_256, PP_NONE, dst, 0, src); }
    void rsqrtps(int dst, RegMem src)
    { vex_nds(0x52, MM_0F, W_IGNORE, L_256, PP_NONE, dst, 0, src); }
    void maxps(int dst, int src0, RegMem src1)
    { vex_nds(0x5f, MM_0F, W_IGNORE, L_256, PP_NONE, dst, src0, src1); }
    void minps(int dst, int src0, RegMem src1)
    { vex_nds(0x5d, MM_0F, W_IGNORE, L_256, PP_NONE, dst, src0, src1); }
    void blendvps(int dst, int src0, RegMem src1, int mask) // dst = mask_hi_bit ? src0 : src1
    { vex_nds(0x4a, MM_0F3A, W_0, L_256, PP_66, dst, src0, src1); emit(mask << 4); }
    void roundps(int dst, RegMem src, round mode)
    { vex_nds(0x08, MM_0F3A, W_IGNORE, L_256, PP_66, dst, 0, src); emit((int)mode); }
    void cvttps2dq(int dst, RegMem src)
    { vex_nds(0x5b, MM_0F, W_IGNORE, L_256, PP_F3, dst, 0, src); }
    void cvtdq2ps(int dst, RegMem src)
    { vex_nds(0x5b, MM_0F, W_IGNORE, L_256, PP_NONE, dst, 0, src); }
    void cmpps(int dst, int src0, RegMem src1, test mode)
    { vex_nds(0xc2, MM_0F, W_IGNORE, L_256, PP_NONE, dst, src0, src1); emit((int)mode); }
    void andps(int dst, int src0, RegMem src1)
    { vex_nds(0x54, MM_0F, W_IGNORE, L_256, PP_NONE, dst, src0, src1); }
    void andnps(int dst, int src0, RegMem src1) // dst = src1 & ~src0
    { vex_nds(0x55, MM_0F, W_IGNORE, L_256, PP_NONE, dst, src0, src1); }
    void orps(int dst, int src0, RegMem src1)
    { vex_nds(0x56, MM_0F, W_IGNORE, L_256, PP_NONE, dst, src0, src1); }
    void xorps(int dst, int src0, RegMem src1)
    { vex_nds(0x57, MM_0F, W_IGNORE, L_256, PP_NONE, dst, src0, src1); }
    void testps(int src0, RegMem src1)
    { vex_nds(0x0e, MM_0F38, W_0, L_256, PP_66, src0, 0, src1); }
    void insertf128(int dst, int into, RegMem from, int hi)
    { vex_nds(0x18, MM_0F3A, W_0, L_256, PP_66, dst, into, from); emit(!!hi); }
    void extractf128(RegMem dst, int src, int hi)
    { vex_nds(0x19, MM_0F3A, W_0, L_256, PP_66, src, 0, dst); emit(!!hi); }
    void broadcastss(int dst, RegMem src)
    { vex_nds(0x18, MM_0F38, W_0, L_256, PP_66, dst, 0, src); }
    void maskmovps(int dst, int mask, RegMem src)
    { vex_nds(0x2c, MM_0F38, W_0, L_256, PP_66, dst, mask, src); }
    void maskmovps_store(RegMem dst, int mask, int src)
    { vex_nds(0x2e, MM_0F38, W_0, L_256, PP_66, src, mask, dst); }
    void movaps(int dst, RegMem src)
    { vex_nds(0x28, MM_0F, W_IGNORE, L_256, PP_NONE, dst, 0, src); }
    void movaps_store(RegMem dst, int src)
    { vex_nds(0x29, MM_0F, W_IGNORE, L_256, PP_NONE, src, 0, dst); }

    // Just a tiny handful of scalar instructions
    void int3() { emit(0xcc); }
    void push(int r) { rex(0, 0, r); emit(0x50 + (r & 7)); }
    void pop(int r) { rex(0, 0, r); emit(0x58 + (r & 7)); }
    void call(RegMem rm)
    { rex(0, rm); emit(0xff); regmem(2, rm); }
    void test32(RegMem a, int imm32) { rex(0, a); emit(0xf7); regmem(0, a); emit32(imm32); }
    void test64(int dst, RegMem src)
    { rexw(dst, src); emit(0x85), regmem(dst, src); }
    void and64_imm8(int reg, char imm)
    { rexw(reg, 0, 0); emit(0x83); regmem(4, reg); emit(imm); }
    void add64_imm32(int dst, int imm)
    { rexw(0, 0, dst); emit(0x81); regmem(0, dst); emit32(imm); }
    void sub64_imm32(int dst, int imm)
    { rexw(0, 0, dst); emit(0x81); regmem(5, dst); emit32(imm); }
    void cmp32_imm8(RegMem lhs, int rhs)
    { rex(0, lhs); emit(0x83); regmem(7, lhs); emit(rhs); }
    void lea(int reg, RegMem ptr)
    { rexw(reg, ptr); emit(0x8d); regmem(reg, ptr); }
    void inc(int reg)
    { rexw(0, 0, reg); emit(0xff); regmem(0, reg); }
    void xor64(int dst, RegMem src)
    { rexw(dst, src); emit(0x33), regmem(dst, src); }
    void not64(RegMem rm)
    { rexw(0, rm); emit(0xf7); regmem(2, rm); }
    void or64(int dst, RegMem src)
    { rexw(dst, src); emit(0x0b), regmem(dst, src); }
    void cmp64(int dst, RegMem src)
    { rexw(dst, src); emit(0x3b), regmem(dst, src); }
    void dec64(RegMem rm)
    { rexw(0, rm); emit(0xff); regmem(1, rm); }

    // Loads and stores come in many delicious and non-orthogonal
    // flavors... Rather than introduce yet more confusion, we
    // preserve the Intel naming here where possible: so no MOVZX for
    // a 32 bit load because that is the specified behavior of regular
    // MOV, the 32->64 sign extension is MOVSXD not MOVSX, etc...  The
    // only non-intel naming is the size and _store suffixes required
    // to disambiguate the C++ function names.
    void movsx8(int dst, RegMem src)
    { rexw(dst, src); emit(0x0f); emit(0xbe); regmem(dst, src); }
    void movsx16(int dst, RegMem src)
    { rexw(dst, src); emit(0x0f); emit(0xbf); regmem(dst, src); }
    void movsxd32(int dst, RegMem src)
    { rexw(dst, src); emit(0x63); regmem(dst, src); }
    void movzx8(int dst, RegMem src)
    { rexw(dst, src); emit(0x0f); emit(0xb6); regmem(dst, src); }
    void movzx16(int dst, RegMem src)
    { rexw(dst, src); emit(0x0f); emit(0xb7); regmem(dst, src); }
    void mov32(int dst, RegMem src)
    { rex(dst, src); emit(0x8b); regmem(dst, src); }
    void mov64(int dst, RegMem src)
    { rexw(dst, src); emit(0x8b); regmem(dst, src); }
    void mov64_imm32(int dst, int imm)
    { rexw(0, 0, dst); emit(0xc7); regmem(0, dst); emit32(imm); }
    void mov8_store(RegMem dst, int src)
    { rex(src, dst); emit(0x88); regmem(src, dst); }
    void mov16_store(RegMem dst, int src)
    { emit(0x66); rex(src, dst); emit(0x89); regmem(src, dst); }
    void mov32_store(RegMem dst, int src)
    { rex(src, dst); emit(0x89); regmem(src, dst); }
    void mov64_store(RegMem dst, int src)
    { rexw(src, dst); emit(0x89); regmem(src, dst); }

    // CALL and Jcc instructions (only a subset are defined here)
    // return a fixup offset for forward displacements.  Note that the
    // address they take is absolute relative to our code start; the
    // conversion to a relative disp is automatic.
    //
    // OPTIMIZATION FIXME: these are 32 bit displacement forms.  It
    // would be nice to support the single-byte form for small jumps,
    // but variable length instructions would require a much more
    // elaborate fixup mechanism...
    int jc(int addr)  { return jcc(0x82, addr); }
    int jnc(int addr) { return jcc(0x83, addr); }
    int jz(int addr)  { return jcc(0x84, addr); } // alias: JE
    int jnz(int addr) { return jcc(0x85, addr); } // alias: JNE
    int jns(int addr)  { return jcc(0x89, addr); }
    int jmp(int addr)  { emit(0xe9); emit32(addr-(size()+4)); return size()-4; }
    void ret() { emit(0xc3); }

    void data8(int b) { emit(b); }
    void data32(int d) { emit32(d); }

    // FIXME: need a way to fixup RIP-relative addresses too
    void fixup(int addr, int val);

    RegMem mem() { return RegMem(); }
    RegMem mem(int base) { return RegMem().base(base); }

    // Returns a comma-separated list of features (currently "avx" and
    // "fma") supported in the local execution environment.
    static std::string detect();

private:
    // enumerants for all the funny opcode fields
    enum vex_w { W_0=0, W_1=1, W_IGNORE=-1 };
    enum vex_l { L_128=0, L_256=1, L_IGNORE=0, L_LZ=0 };
    enum vex_pp { PP_NONE=0, PP_66=1, PP_F3=2, PP_F2=3 }; // Implied SIMD opcode prefixes
    enum vex_mm { MM_0F=1, MM_0F38=2, MM_0F3A=3 }; // Implied prefixes in 3-byte VEX instructions

    void emit(int b) { _code.push_back((unsigned char)b); }
    void emit32(int i) { emit(i); emit(i>>8); emit(i>>16); emit(i>>24); }
    int jcc(int code, int addr);
    void vex(int opcode, int r, int x, int b, vex_mm mm, vex_w w, int vv, vex_l l, vex_pp pp);
    void regmem(int regop, int base, int idx, int scale, int disp);
    void regmem(int regop, RegMem rm);
    void vex_nds(int opcode, vex_mm mm, vex_w w, vex_l l, vex_pp pp, int reg, int vvvv, RegMem rm);

    void rex(int w, int r, int x, int b);
    void rex(int r, int x, int b) { rex(0, r, x, b); }
    void rexw(int r, int x, int b) { rex(1, r, x, b); }

    void rex_(int w, int r, RegMem rm) {
        if(rm._reg == REG_NONE) rex(w, r, rm._idx, rm._base);
        else rex(w, r, 0, rm._reg);
    }
    void rex(int r, RegMem rm) { rex_(0, r, rm); }
    void rexw(int r, RegMem rm) { rex_(1, r, rm); }


    std::vector<unsigned char>& _code;
};

#endif // _AVXGEN_HPP
