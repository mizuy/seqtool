BIN = ../env/bin

PASSWORD = `cat mysql_password`

all: table dbtss

table:
	$(BIN)/seqtooldb load --host=localhost --password=$(PASSWORD)
	touch table

table-ucsc:
	$(BIN)/seqtooldb load --host=genome-mysql.cse.ucsc.edu
	touch table

dbtss:
	$(BIN)/seqtooldb bed --bed_dir=./dbtss/adult-tissue/
	touch dbtss

clear:
	$(BIN)/seqtooldb clear

.PHONY:	table dbtss clear table-ucsc
