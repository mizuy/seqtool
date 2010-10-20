all: example/primers.html example/b2m.html example/actb.html example/btg3.html
	open $+

example/primers.html: example/primers.txt
	primers $< -o $@

output/bisearch.html: example/input.fasta
	bisearch -m normal $< -o output/bisearch.seqv
	seqview output/bisearch.seqv -o $@
	open $@

clean:
	rm -f **/*~
	rm -f #*
	rm -f **/*.pyc
	rm example/*.html

test:
	nosetests

.SUFFIXES: .seqv .html

.seqv.html:
	seqview $< -o $@

.PHONY: test

