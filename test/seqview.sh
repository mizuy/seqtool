#/bin/zsh
rm $1.html
seqview $1.seqv -o $1.html
open $1.html
