all: output/primers.html output/b2m.html output/actb.html

output/primers.html: example/primers.txt
	python -m seqtool.primers $< -o $@
	open $<

output/b2m.html: example/b2m.seqv
	python -m seqtool.seqview $< -o $@

output/actb.html: example/actb.seqv
	python -m seqtool.seqview $+ -o $@

output/bisearch.html: example/input.fasta
	python -m bisearch -m normal $< -o output/bisearch.seqview
	python -m seqtool.seqview output/bisearch.seqview -o $@
	open $@

clean:
	rm -f **/*~
	rm -f #*
	rm -f **/*.pyc

test:
	nosetests

.PHONY: test
