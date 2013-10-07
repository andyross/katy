CFLAGS = -g -Wall -Werror -Os

# Turn on more gcc warnings, but the partial struct initializer syntax
# is something the code relies on to get zero-initialization, so turn
# that warning off.
CFLAGS += -Wextra -Wno-missing-field-initializers

# gcov flags:
ifeq ($(GCOV),1)
CFLAGS += -fprofile-arcs -ftest-coverage
LDFLAGS += -fprofile-arcs
endif

all: vgtest txtest mathtest rastertest

test: vgtest txtest rastertest
	./vgtest --opt-both
	./txtest --opt-both
	@for tool in $^; do \
	  for t in `./$$tool --tests`; do \
	    echo Register check \(unoptimized\): $$t; \
	    ./$$tool --no-optimize $$t --dump-asm | ./checkasm.pl $$t; \
	    echo Register check: $$t; \
	    ./$$tool $$t --dump-asm | ./checkasm.pl $$t; \
	  done; \
	done

vgtest: vgtest.o test.o vecgen.o vecgen-avx.o avxgen.o vecgen-opt.o math.o
	g++ $(LDFLAGS) -o $@ $^

txtest: txtest.o test.o texture.o vecgen.o vecgen-avx.o avxgen.o vecgen-opt.o
	g++ $(LDFLAGS) -o $@ $^

rastertest: rastertest.o test.o texture.o vecgen.o vecgen-avx.o avxgen.o vecgen-opt.o raster.o framebuf.o
	g++ $(LDFLAGS) -o $@ $^ -ljpeg

mathtest: mathtest.o vecgen.o vecgen-avx.o avxgen.o vecgen-opt.o math.o
	g++ $(LDFLAGS) -o $@ $^ -pthread -lm

# Generate a quick gcov report of unhit lines
coverage:
	$(MAKE) clean
	$(MAKE) GCOV=1 coverage_internal
	rm -f *.o

coverage_internal: txtest vgtest rastertest
	-./vgtest --log --dump-code --dump-asm --opt-both --no-halt > /dev/null
	-./txtest --log --dump-code --dump-asm --opt-both --no-halt > /dev/null
	-./rastertest --log --dump-code --dump-asm --opt-both --no-halt > /dev/null
	sh -c 'for i in *pp; do gcov $$i; done' > /dev/null 2>&1
	sh -c 'for i in *pp; do echo; echo $$i:; fgrep "#####" $$i.gcov; done; true' 2>/dev/null

clean:
	rm -f vgtest txtest mathtest mpqtest rastertest *.d *.o *.gc?? tri1.ppm

%.o : %.cpp
	g++ -I. -c -MMD $(CFLAGS) -o $@ $<

-include *.d
