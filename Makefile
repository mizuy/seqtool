CPP = clang++
CFLAGS = -Wall -std=c++0x -stdlib=libc++ -L/usr/local/lib -I/usr/local/include

all: seqview tssview geneview get_genbank primers bisearch
	
seqview:
	seqview input/seqview.seqv -o _output/seqview.html

tssview:
	tssview input/tssview.tssv -o _output/tssview.html

geneview:
	geneview TP53 -o _output/geneview.html

get_genbank:
	get_genbank TP53 > _output/TP53.gb

primers:
	primers input/primers.txt -o _output/primers.html 
	primers input/primers.txt -c -o _output/primers.csv 

bisearch: bin/bisearch
	./bin/bisearch input/test.fasta > _output/bisearch.seqv

bin/bisearch: c/bisearch.cpp c/bisearch.h c/nucleotide.h c/main.cpp
	$(CPP) $(CFLAGS) -O1 -lboost_program_options c/bisearch.cpp c/main.cpp -o bin/bisearch


clean:
	rm -f **/*~
	rm -f #*
	rm -f **/*.pyc
	rm example/*.html

test:
	nosetests

.PHONY: test all

