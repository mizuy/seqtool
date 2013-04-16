CPP = clang++
CFLAGS = -Wall -std=c++0x -stdlib=libc++ -L/usr/local/lib -I/usr/local/include

all: build

build:
	cd input; make build

bin/bisearch: c/bisearch.cpp c/bisearch.h c/nucleotide.h c/main.cpp
	$(CPP) $(CFLAGS) -O1 -lboost_program_options c/bisearch.cpp c/main.cpp -o bin/bisearch

clean:
	rm -f **/*~
	rm -f #*
	rm -f **/*.pyc
	rm example/*.html

test:
	bin/nosetests

pdf: source.pdf
source.pdf: seqtool/*.py seqtool/**/*.py
	enscript seqtool/*.py seqtool/**/*.py --font=Courier6 --highlight=python --line-numbers --landscape --columns=2 --color -o source.ps
	pstopdf source.ps -o source.pdf
	rm source.ps

.PHONY: test all build clean

