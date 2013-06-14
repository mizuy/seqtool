
def remove_space(d):
    ret = {}
    for k,v in list(d.items()):
        ret[k] = v.replace(' ','')
    return ret

'''
primer=unambiguous or ambiguous(any), template=unambiguous
 basematch_unambiguous

primer=unambiguous, template=ambiguous(any)

primer=ambiguous(any), template=ambiguous(any)
 basematch_partial

primer=ambiguous(all), template=ambiguous(any)
'''


basematch_unambiguous = remove_space({
     'A': 'A   ',
     'T': ' T  ',
     'G': '  G ',
     'C': '   C',
     'W': 'AT  ',
     'R': 'A G ',
     'M': 'A  C',
     'K': ' TG ',
     'Y': ' T C',
     'S': '  GC',
     'B': ' TGC',
     'V': 'A GC',
     'H': 'AT C',
     'D': 'ATG ',
     'N': 'ATGC'})

basematch_partial = remove_space({
         # ATGC WRMKYS BVHD N
     'A': 'A    WRM     VHD N',
     'T': ' T   W  KY  B HD N',
     'G': '  G   R K S BV D N',
     'C': '   C   M YS BVH  N',
     'W': 'AT   WRMKY  BVHD N',
     'R': 'A G  WRMK S BVHD N',
     'M': 'A  C WRM YS BVHD N',
     'K': ' TG  WR KYS BVHD N',
     'Y': ' T C W MKYS BVHD N',
     'S': '  GC  RMKYS BVHD N',
     'B': ' TGC WRMKYS BVHD N',
     'V': 'A GC WRMKYS BVHD N',
     'H': 'AT C WRMKYS BVHD N',
     'D': 'ATG  WRMKYS BVHD N',
     'N': 'ATGC WRMKYS BVHD N'})

basematch_subset = remove_space({
         # ATGC WRMKYS BVHD N
     'A': 'A    WRM     VHD N',
     'T': ' T   W  KY  B HD N',
     'G': '  G   R K S BV D N',
     'C': '   C   M YS BVH  N',
     'W': '     W        HD N',
     'R': '      R      V D N',
     'M': '       M     VH  N',
     'K': '        K   B  D N',
     'Y': '         Y  B H  N',
     'S': '          S BV   N',
     'B': '            B    N',
     'V': '             V   N',
     'H': '              H  N',
     'D': '               D N',
     'N': '                 N'})

basematch_superset = remove_space({
         # ATGC WRMKYS BVHD N
     'A': 'A                 ',
     'T': ' T                ',
     'G': '  G               ',
     'C': '   C              ',
     'W': 'AT   W            ',
     'R': 'A G   R           ',
     'M': 'A  C   M          ',
     'K': ' TG     K         ',
     'Y': ' T C     Y        ',
     'S': '  GC      S       ',
     'B': ' TGC    KYS B     ',
     'V': 'A GC  RM  S  V    ',
     'H': 'AT C W M Y    H   ',
     'D': 'ATG  WR K      D  ',
     'N': 'ATGC WRMKYS BVHD N'})

def oligo_regex(seq, match=basematch_partial):
    return ''.join(['[{}]'.format(match[s]) for s in str(seq).upper()])

def base_match(base_i, base_j, match=basematch_partial):
    return base_j in match[base_i]
