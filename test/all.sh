seqview input/seqview.seqv -o _output/seqview.html
tssview input/tssview.tssv -o _output/tssview.html
geneview TP53 -o _output/geneview.html
get_genbank TP53 > _output/TP53.gb
primers input/primers.txt -o _output/primers.html 
primers input/primers.txt -c -o _output/primers.csv 
