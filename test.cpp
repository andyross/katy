// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include <sys/mman.h>
#include <cstdlib>
#include <set>
#include "avxgen.hpp"
#include "test.hpp"

using namespace std;

// Booleans settable via command line
bool run_interpreter = true;
bool run_avx = true;
bool random_fill = true;
bool log_execution = false;
bool optimize = true;
bool opt_both = false;
bool dump_asm = false;
bool dump_code = false;
bool halt = true;
bool breakpoint = false;

set<string> tests_enabled;
bool avxmode;
int simd_index;

struct optrec { const char *name; bool *ptr; } boolopts[] = {
    { "interpreter", &run_interpreter },
    { "avx", &run_avx },
    { "log", &log_execution },
    { "optimize", &optimize },
    { "opt-both", &opt_both },
    { "dump-asm", &dump_asm },
    { "dump-code", &dump_code },
    { "random-fill", &random_fill },
    { "halt", &halt },
    { "break", &breakpoint },
    {}
};

const unsigned int MAX_TESTS = 1024;
static const struct test_rec *tests[MAX_TESTS];
static unsigned int ntests = 0;

test_rec::test_rec(const char *n, bool(*f)())
    : name(n), fn(f)
{
    if(ntests <= MAX_TESTS)
        tests[ntests++] = this;
    else
        fprintf(stderr, "TOO MANY TESTS, skipping %s\n", n);
}

bool parse_command_line(int argc, char **argv)
{
    for(int i=1; i<argc; i++) {
        string arg = argv[i];
        bool val = true;

        if(arg.substr(0, 2) != "--") {
            bool ok = false;
            for(unsigned int j=0; j<ntests; j++)
                if(arg == tests[j]->name)
                    ok = true;
            if(!ok) {
                printf("No such test: %s\n", argv[i]);
                return false;
            }
            tests_enabled.insert(argv[i]);
            continue;
        }
        arg = arg.substr(2);
        if(arg.substr(0, 3) == "no-") {
            val = false;
            arg = arg.substr(3);
        }

        if(arg == "simd-index" && i < (argc-1) && val) {
            simd_index = atoi(argv[++i]);
        } else if(arg == "tests") {
            // Dump list of tests
            for(unsigned int i=0; i<ntests; i++)
                printf("%s\n", tests[i]->name);
            exit(0);
        } else {
            int j;
            for(j=0; boolopts[j].name; j++) {
                if(arg == boolopts[j].name) {
                    *boolopts[j].ptr = val;
                    break;
                }
            }
            if(!boolopts[j].name)
                return false;
        }
    }
    return true;
}

// The executable memory used by the AVX test
union funcrec {
    void *mem;
    void (*code)(float* consts, float* ins, float* outs, void **mems, float *msk);
} func;

// Allocates a properly aligned SIMD buffer of n registers and fills
// it with either zeros or random garbage
float* simd_mem(int n)
{
    const size_t ymmsz = 8*4;
    float *result = 0;
    if(posix_memalign((void**)&result, ymmsz, n*ymmsz) == 0)
        for(unsigned int i=0; i<(n*ymmsz/sizeof(int)); i++)
            ((int*)result)[i] = random_fill ? rand() : 0;
    return result;
}

string optstring()
{
    string opts;
    if(!optimize) opts += opts.size() ? ",no-optimize" : "no-optimize";
    if(dump_asm)  opts += opts.size() ? ",dump-asm" : "dump-asm";
    if(breakpoint) opts += opts.size() ? ",break" : "break";
    return opts;
}

void vecgen_run(vecgen &vg, float *consts, float *in, float *out, void **mems)
{
    vg.codegen();
    vecgen_run_no_gen(vg, consts, in, out, mems);
}

void vecgen_run_no_gen(vecgen &vg, float *consts, float *in, float *out, void **mems)
{
    if(dump_code)
        vg.dump_code(stdout);

    if(!avxmode) {
        vg.interpret(consts, in, out, mems, log_execution ? stdout : NULL);
        return;
    }

    float *avxin = simd_mem(vg.ninputs());
    float *avxout = simd_mem(vg.noutputs());
    float *msk = simd_mem(1);

    for(int i=0; i<vg.ninputs(); i++)
        avxin[8*i+simd_index] = in[i];
    for(int i=0; i<vg.noutputs(); i++)
        avxout[8*i+simd_index] = out[i]; // Outputs can be left unchanged, so they're "inputs" too

    // MSK gets 0x80000000 in the enabled slot, otherwise zeros
    for(int i=0; i<8; i++)
        ((int*)msk)[i] = (i == simd_index) << 31;

    if(dump_asm)
        printf("Calling generated code at %p\n", func.mem);

    vg.get_code((char*)func.mem);
    func.code(consts, avxin, avxout, mems, msk);

    // Copy back the outputs
    for(int i=0; i<vg.noutputs(); i++)
        out[i] = avxout[8*i+simd_index];

    free(avxin);
    free(avxout);
    free(msk);
}

bool run_all_tests()
{
    for(unsigned int i=0; i<ntests; i++) {
        if(tests_enabled.size())
            if(tests_enabled.find(tests[i]->name) == tests_enabled.end())
                continue;

        int iterations = 0;
        for(int j=0; j<2; j++) {
            avxmode = (j == 1);
            if((avxmode && !run_avx) || (!avxmode && !run_interpreter))
                continue;

            // Don't dump the assembly multiple times
            if(iterations++)
                dump_asm = dump_code = false;

            printf("## TEST: %s (%s%s)\n", tests[i]->name, avxmode ? "AVX" : "interpreter",
                   optimize ? ", opt" : "");
            bool result = tests[i]->fn();
            printf("## %s %s\n", tests[i]->name, result ? "SUCCEEDED" : "FAILED");
            if(!result && halt)
                return false;
        }
    }
    return true;
}

int main(int argc, char** argv)
{
    run_avx = avxgen::detect().find("avx") != string::npos;

    if(!parse_command_line(argc, argv)) {
        printf("Usage: %s [--test <test1_name>] [--test <test2_name>]...\n",
               argv[0]);
        printf("       [--tests]  (prints list of tests to stdout)\n");
        printf("       [--simd-index <0-7>] (default=%d)\n", simd_index);
        for(int i=0; boolopts[i].name; i++)
            printf("       [--{no-}%s] (default=%s)\n",
                   boolopts[i].name,
                   *boolopts[i].ptr ? "true" : "false");
        printf("Available tests:\n");
        for(unsigned int i=0; i<ntests; i++)
            printf("  %s\n", tests[i]->name);
        return 1;
    }

    // FIXME: fixed size of 64k with no overrun detection...
    func.mem = mmap(0, 65536, PROT_EXEC|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
    if(func.mem == MAP_FAILED) {
        fprintf(stderr, "mmap of code region failed\n");
        return 1;
    }

    if(opt_both) {
        optimize = false;
        if(!run_all_tests()) return 1;
        optimize = true;
        if(!run_all_tests()) return 1;
    } else {
        return run_all_tests() ? 0 : 1;
    }

    return 0;
}
