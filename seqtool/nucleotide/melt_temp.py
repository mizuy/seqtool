from math import log10, log, exp

__all__ = ['melting_temperature_unambiguous','melting_temperature_unmethyl']

# gas constant, [J/molK]
R = 1.987
K = 273.15

class PCRMixture(object):
    def __init__(self, c_na, c_k, c_mg, c_primer, c_ntp=None):
        """
        c_na, c_mg: final concentration in molar of Na and Mg in PCR reaction mix
                    these are used for calculation of salt concentration
        c_primer:   final concentration in molar of each primer in PCR reaction mix
        c_ntp:      
        """
        self.c_na = c_na
        self.c_k = c_k
        self.c_mg = c_mg
        self.c_primer = c_primer
        self.c_ntp = c_ntp

    def cation_conc(self):
        """
        Ahsen, von, N.
        Oligonucleotide Melting Temperatures under PCR Conditions: Nearest-Neighbor Corrections for Mg2+, Deoxynucleotide Triphosphate, and Dimethyl Sulfoxide Concentrations with Comparison to Alternative Empirical Formulas.

        c_na + c_k + 120 * (c_mg-c_ntp)**0.5
        """
        return self.c_na + self.c_k + 4*(self.c_mg**0.5)

DEFAULT_C_NA = 33.*10**-3
DEFAULT_C_K = 0
DEFAULT_C_MG = (2+4.5)*10**-3
DEFAULT_C_PRIMER = 0.5*10**-6
DEFAULT_C_NTP = DEFAULT_C_PRIMER

DEFAULT_MIX = PCRMixture(DEFAULT_C_NA, DEFAULT_C_K,  DEFAULT_C_MG, DEFAULT_C_PRIMER, DEFAULT_C_NTP)

    
def melting_temperature_unmethyl(seq, pcr_mix=DEFAULT_MIX, unmethyl=True):
    seq = str(seq).upper()
    if unmethyl:
        seq = seq.replace('R','A').replace('Y','T')
    else:
        seq = seq.replace('R','G').replace('Y','C')

    return melting_temperature_unambiguous(seq, pcr_mix)


def melting_temperature_unambiguous(seq, pcr_mix=DEFAULT_MIX):
    '''
    unmethyl:   if 'unmethyl' is true, all CpGs of template are assumed to be unmethyled
                then unmethyl version of primer are used for calculation
    '''

    return wetmur_tm(str(seq), pcr_mix.cation_conc(), pcr_mix.c_primer)

def complex_fraction(seq, anneal_temp, pcr_mix=DEFAULT_MIX):
    return sl_fraction(str(seq), anneal_temp, pcr_mix.cation_conc(), pcr_mix.c_primer)

##################################################################
## Wetmur, J. G., & Fresco, J. (1991). 
## DNA Probes: Applications of the Principles of Nucleic Acid Hybridization.
##################################################################

"""
Tm = T0 DeltaH / DeltaH - DeltaG + RT ln(Cp) + salt_term - K

salt_term = + 16.6 log( Csalt / 1+0.7*CSalt ) + 3.85
T0 = 298.2 // = 25 + K, room temperature.

where, for Csalt = 1, salt_term = 0
"""

NNT_DH = {
    'AA': 9.1,
    'TT': 9.1,
    'AT': 8.6,
    'TA': 6.0,
    'CA': 5.8,
    'TG': 5.8,
    'GT': 6.5,
    'AC': 6.5,
    'CT': 7.8,
    'AG': 7.8,
    'GA': 5.6,
    'TC': 5.6,
    'CG': 11.9,
    'GC': 11.1,
    'GG': 11.0,
    'CC': 11.0,
}
NNT_DG = {
    'AA': 1.55,
    'TT': 1.55,
    'AT': 1.25,
    'TA': 0.85,
    'CA': 1.15,
    'TG': 1.15,
    'GT': 1.40,
    'AC': 1.40,
    'CT': 1.45,
    'AG': 1.45,
    'GA': 1.15,
    'TC': 1.15,
    'CG': 3.05,
    'GC': 2.70,
    'GG': 2.3,
    'CC': 2.3,
}


def wetmur_tm(seqstr, cation_conc, c_primer):
    l = len(seqstr)
    t0 = 298.2 # room temp.
    d_h_e = 5.
    d_h_p = -1000. * ( 2*d_h_e + sum(NNT_DH[seqstr[i:i+2]] for i in range(l-1)) )
    d_g_e = 1.
    d_g_i = -2.2
    d_g_p = -1000. * ( 2*d_g_e + d_g_i + sum(NNT_DG[seqstr[i:i+2]] for i in range(l-1)) )
    salt_term =  16.6*log10(cation_conc / (1.+0.7*cation_conc)) + 3.85
    t_p = t0*d_h_p / (d_h_p-d_g_p + R*t0*log(c_primer) ) + salt_term - K
    return t_p
    
    

##################################################################
## Xia, T., SantaLucia, J., Burkard, M. E., Kierzek, R., Schroeder, S. J., Jiao, X., Cox, C., et al. (1998).
## Thermodynamic parameters for an expanded nearest-neighbor model for formation of RNA duplexes with Watson-Crick base pairs.
##################################################################

def sl_delta_s(seqstr, cation_conc):
    l = len(seqstr)
    ret = sum([SantaLuciaS[seqstr[i:i+2]] for i in range(l-1)])
    # terminal modification
    ret += SantaLuciaTerminalS[seqstr[0]]
    ret += SantaLuciaTerminalS[seqstr[-1]]
    delta_s = -1 * ret

    """
    SantaLucia equation doesnt have salt_term but have modified DeltaS for caion_conc adjustment.
    """
    return delta_s - 0.368 * (len(seqstr)-1) * log10(cation_conc)

def sl_delta_h(seqstr):
    l = len(seqstr)
    ret = sum([SantaLuciaH[seqstr[i:i+2]] for i in range(l-1)])
    # terminal modification
    ret += SantaLuciaTerminalH[seqstr[0]]
    ret += SantaLuciaTerminalH[seqstr[-1]]
    delta_h = -1 * ret

    return delta_h

def sl_delta_g(seqstr, cation_conc, temp_k):
    """
    DeltaG = 1000 * DeltaH - T * DeltaS
    """
    return 1000. * sl_delta_h(seqstr) - temp_k * sl_delta_s(seqstr, cation_conc)

def sl_tm(seqstr, cation_conc, primer_conc):
    """
    Tm = 1000. * DeltaH / [DeltaS - R log(primer_conc / 4)]
    """
    delta_s = sl_delta_s(seqstr, cation_conc) 
    delta_h = sl_delta_h(seqstr)
    tm = 1000. * delta_h / ( delta_s - R * log10(primer_conc/4.) )
    return tm - K

def sl_fraction(seqstr, aneal_temp, cation_conc, primer_conc):
    """
    A novel strategy to design highly specific PCR primers based on the stability and uniqueness of 3'-end subsequences.
    Miura, F. (2005). Bioinformatics, 21(24), 4363-4370.

    fraction = conc of primer-template complex / conc of template
    CK = primer_conc * Kas
    Kas = exp(deltaG / RT)
    fraction = CK / ( 1+CK )
    """
    temp_k = aneal_temp + K
    delta_g = sl_delta_g(seqstr, cation_conc, temp_k)
    a = primer_conc * exp(delta_g / (R * temp_k))
    f = a/(1.+a)
    return f


SantaLuciaTerminalS = {
'A' : 4.1,
'T' : 4.1,
'C' : -2.8,
'G' : -2.8,
}

SantaLuciaTerminalH = {
'A' : 2.3,
'T' : 2.3,
'C' : 0.1,
'G' : 0.1,
}
    
SantaLuciaH = {
'AA' :  -7.9,
'AC' :  -8.4,
'AG' :  -7.8,
'AT' :  -7.2,
'CA' :  -8.5,
'CC' :  -8.0,
'CG' :  -10.6,
'CT' :  -7.8,
'GA' :  -8.2,
'GC' :  -9.8,
'GG' :  -8.0,
'GT' :  -8.4,
'TA' :  -7.2,
'TC' :  -8.2,
'TG' :  -8.5,
'TT' :  -7.9,
}

SantaLuciaS = {
'AA' :  -22.2,
'AC' :  -22.4,
'AG' :  -21.0,
'AT' :  -20.4,
'CA' :  -22.7,
'CC' :  -19.9,
'CG' :  -27.2,
'CT' :  -21.0,
'GA' :  -22.2,
'GC' :  -24.4,
'GG' :  -19.9,
'GT' :  -22.4,
'TA' :  -21.3,
'TC' :  -22.2,
'TG' :  -22.7,
'TT' :  -22.2,
}
