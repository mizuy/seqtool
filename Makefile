CPP = clang++
CFLAGS = -Wall -std=c++0x -stdlib=libc++ -L/usr/local/lib -I/usr/local/include

all: seqview tssview geneview get_genbank primers bin/bisearch

build: all

sequencing:
#	bin/sequencing input/foxa2 -t input/foxa2.fasta
#	bin/sequencing input/oxbs -t input/oxbscontrol.fasta

#	bin/sequencing input/seq_test -t input/oxbscontrol.fasta

seqview:
	bin/seqview input/seqview.seqv -o _output/seqview.html

tssview:
	bin/tssview input/tssview.tssv -o _output/tssview.html

geneview:
	bin/geneview TP53 -o _output/geneview.html

get_genbank:
	bin/get_genbank TP53 > _output/TP53.gb

primers:
	bin/primers input/primers.txt -o _output/primers.html 
	bin/primers input/primers.txt -c -o _output/primers.csv 

bisearch: bin/bisearch
	bin/bisearch input/test.fasta > _output/bisearch.seqv
	bin/seqview _output/bisearch.seqv -o _output/bisearch.html

bin/bisearch: c/bisearch.cpp c/bisearch.h c/nucleotide.h c/main.cpp
	$(CPP) $(CFLAGS) -O1 -lboost_program_options c/bisearch.cpp c/main.cpp -o bin/bisearch

clean:
	rm -f **/*~
	rm -f #*
	rm -f **/*.pyc
	rm example/*.html

test:
	bin/nosetests

pdf:
	enscript seqtool/*.py seqtool/**/*.py --font=Courier6 --highlight=python --line-numbers --landscape --columns=2 --color -o source.ps

.PHONY: test all

