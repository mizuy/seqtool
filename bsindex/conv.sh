set -x verbose

for f in chromFa/*.fa
do
    echo "processing $f file..."
    ./convert p m $f > bsfa/`basename $f .fa`_pm.fa
    ./convert p u $f > bsfa/`basename $f .fa`_pu.fa
    ./convert n m $f > bsfa/`basename $f .fa`_nm.fa
    ./convert n u $f > bsfa/`basename $f .fa`_nu.fa
done
