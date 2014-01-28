SHELL = /bin/sh
CFLAGS = -O2 -DNDEBUG
VERSION = `grep 'define BEAND_VERSION_STRING' beandm.cpp | sed 's/.*"\(.*\)"/\1/'`

bns: beandm.cpp lognum.hpp
	g++ $(CFLAGS) beandm.cpp -o beandmulti -lgsl -lgslcblas -lm -lboost_program_options 

clean:
	rm -f beandmulti beandm.tar.gz gmon.out

recompile: clean beandm

package: clean
	tar -cpzf BEANDiscoMulti-$(VERSION).tar.gz *

