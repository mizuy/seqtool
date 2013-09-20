BIN = env/bin

all: build

dbload:
	cd database; make

bisearch:
	cd c; make

clean:
	find . -name '*~' -delete
	find . -name '*.pyc' -delete
	find . -name '__pycache__' -delete

cleanall: clean
	rm -rf develop-egg parts *.egg-info dist

coverage:
	$(BIN)/nosetests --with-coverage --cover-html --with-doctest --cover-package=seqtool

test:
	$(BIN)/nosetests --with-doctest

ctags:
	ctags -e -R seqtool

examples:
	cd example; make

build: bisearch
	cd input; make build

bootstrap:
	virtualenv --python=python3.3 --system-site-packages env
	python setup.py develop

pdf: source.pdf
source.pdf: seqtool/*.py seqtool/**/*.py
	enscript seqtool/*.py seqtool/**/*.py --font=Courier6 --highlight=python --line-numbers --landscape --columns=2 --color -o source.ps
	pstopdf source.ps -o source.pdf
	rm source.ps

.PHONY: test coverage all build clean examples cleanall dbload bsearch bootstrap

