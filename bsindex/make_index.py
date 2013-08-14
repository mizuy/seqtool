import glob

build = 'time ~/work/bowtie-1.0.0/bowtie-build '

def a():
    print build + ','.join(glob.glob('bsfa/*_pm.fa')) + ' bs_pos_met'
    print build + ','.join(glob.glob('bsfa/*_pu.fa')) + ' bs_pos_unm'
    print build + ','.join(glob.glob('bsfa/*_nm.fa')) + ' bs_neg_met'
    print build + ','.join(glob.glob('bsfa/*_nu.fa')) + ' bs_neg_unm'

    
a()
