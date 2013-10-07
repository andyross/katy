// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include <stdexcept>
#include "avxgen.hpp"

using namespace std;

// r - Top/REX bit of ModRM reg field
// x - Top/REX bit of SIB index field
// b - Top/REX bit of ModRM r/m or SIB base field
// mm - Leading opcode bytes
// w - opcode-specific extra bit
// vv - non-destructive (NDS/NDD) register specifier, zero if unused
// l - vector length (1 == 256 bit)
// pp - implied SIMD prefix
void avxgen::vex(int opcode, int r, int x, int b, vex_mm mm, vex_w w, int vv, vex_l l, vex_pp pp)
{
    // FIXME: docs are a little ambiguous about when the 2-byte form
    // can be used.  The mmmmm field definitely defaults to the "0F"
    // form, and I *think* the X and B bits default to "bottom half of
    // the register set" (i.e. as if there were no REX prefix), but
    // I'm assuming W must be "WIG" or else we need an explicit field
    // because I can't find anything saying that it defaults to 1 or
    // 0.
    if(w == W_IGNORE && mm == MM_0F && !x && !b) {
        // 2-byte VEX
        emit(0xc5);
        emit((!r << 7) | ((~vv & 0xf) << 3) | (l << 2) | (pp & 3));
    } else {
        // 3-byte VEX
        // (note: need separate WIG tag to be able to choose 2-byte, but must
        // encode w field as 0 when using 3-byte)
        if(w == W_IGNORE)
            w = W_0;
        emit(0xc4);
        emit((!r << 7) | (!x << 6) | (!b << 5) | (mm & 0x1f));
        emit((!!w << 7) | ((~vv & 0xf) << 3) | (l << 2) | (pp & 3));
    }
    emit(opcode);
}

// Computes and emits ModRM/SIB bytes to produce
// [base+idx<<scale+disp] addressing.
//
// The scale argument should be one of [0-3].  Pass REG_NONE for base
// and index to indicate none/unused.  Pass REG_RIP for base to
// indicate RIP-relative addressing (in which case idx/scale must be
// unused and disp is an offset relative to the next instruction).
//
// Note: pass in the 4-bit long mode register IDs, even though this
// only encodes the low 3 bits.  The high bit is sometimes significant
// (to distinguish RSP from R12 when used as idx) for encoding
// decisions.
void avxgen::regmem(int regop, int base, int idx, int scale, int disp)
{
    // Correct attempted use of RSP as index and needless SIB requests
    if((idx == 4 && !scale) || (base == REG_NONE && !scale))
        swap(base, idx);

    if(idx == 4 || (base == REG_RIP && idx != REG_NONE))
        throw domain_error("unencodable ModRM/SIB values");

    // Figure out the basic stuff
    bool sib = idx != REG_NONE || base == REG_NONE;
    int mod = disp ? ((disp < -128 || disp > 127) ? 2 : 1) : 0;
    int db = mod == 2 ? 4 : mod; // disp bytes

    // Some situations need special encodings, map them away from their "natural" form:
    if((base & 7) == 4) {
        sib = true; // base==4 means "use SIB", not "RSP/R12"
        if(idx == REG_NONE)
            idx = 4;
    } else if((base & 7) == 5 && mod == 0)
        mod = 1; // Can't encode [RBP] or [R13] with no disp.  Use a disp8.

    // ... because other special forms get mapped into those slots:
    if(base == REG_RIP) {
        mod = 0; base = 5; db = 4; sib = false;
    } else if(base == REG_NONE) {
        mod = 0; base = 5; db = 4; sib = true;
        if(idx == REG_NONE)
            idx = 4;
    }

    emit((mod << 6) | ((regop & 7) << 3) | ((sib ? 4 : base) & 7));
    if(sib)
        emit((scale << 6) | ((idx & 7) << 3) | (base & 7));
    for(int i=0; i<db; i++)
        emit((disp >> (8*i)) & 0xff);
}

void avxgen::regmem(int regop, RegMem rm)
{
    if(rm._reg == REG_NONE)
        regmem(regop, rm._base, rm._idx, rm._scale, rm._disp);
    else
        emit((3 << 6) | ((regop & 7) << 3) | (rm._reg & 7)); // ModRM
}


void avxgen::vex_nds(int opcode, vex_mm mm, vex_w w, vex_l l, vex_pp pp, int reg, int vvvv, RegMem rm)
{
    if(rm._reg == REG_NONE) {
        int idx = rm._idx == REG_NONE ? 0 : rm._idx;
        vex(opcode, !!(reg & 8), !!(idx & 8), !!(rm._base & 8), mm, w, vvvv, l, pp);
        regmem(reg, rm);
    } else {
        vex(opcode, !!(reg & 8), 0, !!(rm._reg & 8), mm, w, vvvv, l, pp);
        emit((3 << 6) | ((reg & 7) << 3) | (rm._reg & 7)); // ModRM
    }
}

void avxgen::fixup(int addr, int val)
{
    _code[addr]   = (unsigned char)val;
    _code[addr+1] = (unsigned char)(val>>8);
    _code[addr+2] = (unsigned char)(val>>16);
    _code[addr+3] = (unsigned char)(val>>24);
}

int avxgen::jcc(int code, int addr)
{
    int next = size() + 6;
    emit(0x0f);
    emit(code);
    emit32(addr-next);
    return size() - 4;
}

void avxgen::rex(int w, int r, int x, int b)
{
    r = r == REG_NONE ? 0 : !!(r & 8);
    x = x == REG_NONE ? 0 : !!(x & 8);
    b = b == REG_NONE ? 0 : !!(b & 8);
    if(w | r | x | b)
        emit(0x40 | ((!!w)<<3) | (r<<2) | (x<<1) | b);
}

#if defined(__x86_64) || defined(__i686)
static void cpuid(int key, unsigned int *a, unsigned int *b, unsigned int *c, unsigned int *d)
{ asm("cpuid" : "=a"(*a), "=b"(*b), "=c"(*c),"=d"(*d) : "0"(key)); }
static void xgetbv(int key, unsigned int *eax, unsigned int *edx)
{ asm("xgetbv" : "=a"(*eax), "=d"(*edx) : "c"(key)); }
#else
static void cpuid(int key, unsigned int *a, unsigned int *b, unsigned int *c, unsigned int *d)
{ *a = *b = *c = *d = 0; }
static void xgetbv(int key, unsigned int *eax, unsigned int *edx)
{ *eax = *edx = 0; }
#endif

string avxgen::detect()
{
    string s = "avx";
    unsigned int a, b, c, d;

    // AVX bit in CPUID:01H
    cpuid(1, &a, &b, &c, &d);
    if(!(c & (1<<28)))
        return "";

    // FMA bit in the same CPUID:01H word
    if(c & (1<<12))
        s += ",fma";

    // OXSAVE bit: indicates support for XGETBV instruction ...
    if(!(c & (1<<27)))
        return "";

    // ... which returns OS support for the extended registers
    xgetbv(0, &a, &d);
    if(!(a & (1<<1)))
        return ""; // No XMM/SSE support in OS
    if(!(a & (1<<2)))
        return ""; // No YMM/AVX support in OS

    return s;
}
