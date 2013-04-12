
def dna_mw(length, ds=True):
    mw = 330 * length
    if ds:
        mw = 2 * mw
    return mw

def dna_ug_to_pmol(weight, length, ds=True):
    return 1.* 10**6 * weight / dna_mw(length, ds)

def dna_pmol_to_ug(mol, length, ds=True):
    return 1.* mol * dna_mw(length, ds) / 10**6

print dna_pmol_to_ug(dna_ug_to_pmol(50, 3000), 3000), 50
print 'vector', dna_ug_to_pmol(50, 3000), 'pmol'
print 'insert', dna_pmol_to_ug(dna_ug_to_pmol(50, 3000) * 0.1, 338), 'ug'
