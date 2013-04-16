import glob

build = 'time bowtie-build '

def a():
    print build + ','.join(glob.glob('*_bs_sense_met.fa')) + ' bs_sense_met'
    print build + ','.join(glob.glob('*_bs_sense_unmet.fa')) + ' bs_sense_unmet'
    print build + ','.join(glob.glob('*_bs_antisense_met.fa')) + ' bs_antisense_met'
    print build + ','.join(glob.glob('*_bs_antisense_unmet.fa')) + ' bs_antisense_unmet'

    
"nohup ./build.sh > out.log 2> err.log < /dev/null &"
