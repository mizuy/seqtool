CPP = clang++
CFLAGS = -Wall -std=c++0x -stdlib=libc++ -L/usr/local/lib -I/usr/local/include

all: build

build:
	cd input; make build

bin/bisearch: c/bisearch.cpp c/bisearch.h c/nucleotide.h c/main.cpp c/primer.cpp
	$(CPP) $(CFLAGS) -O2 -lboost_program_options c/bisearch.cpp c/main.cpp c/primer.cpp -o bin/bisearch

examples:
	bin/seqview example/actb.seqv
	bin/seqview example/b2m.seqv
	bin/seqview example/btg3.seqv


clean:
	rm -f **/*~
	rm -f #*
	rm -f **/*.pyc
	rm example/*.html

cleanall: clean
	rm -rf develop-egg parts *.egg-info dist

test:
	bin/nosetests

dbload:
	cd database; make
	bin/database_load --chrom_tab_file=database/chrom.tab --ucsc_tab_file=database/ucsc.tab --hgnc_tab_file=database/hgnc.tab

pdf: source.pdf
source.pdf: seqtool/*.py seqtool/**/*.py
	enscript seqtool/*.py seqtool/**/*.py --font=Courier6 --highlight=python --line-numbers --landscape --columns=2 --color -o source.ps
	pstopdf source.ps -o source.pdf
	rm source.ps

.PHONY: test all build clean examples cleanall dbload

