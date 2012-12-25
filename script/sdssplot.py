import seqtool.pcr as pcr
import matplotlib.pyplot as pp

l = [x.split(":") for x in """
IRX5-MA:  CGAACGACGCCCTACGAAAACG
IRX5-MS:  CGTCGCGGGTCGGAGTTTC
IRX5-UA:  CAACAAACCAAACAACACCCTACAAA
IRX5-US:  GGAGGTGTTGTGGGTTGGAGTTTTG
IRX5-MS2: CGCGGCGGTCGTAGTC
IRX5-MA3: GACGACCGCGAAATCGACG
IRX5-US2: GAGTGGGTGTGGTGGTTGTAGTTG
IRX5-UA3: AAAAACAACAACAACAACAACCACAA
""".strip().split("\n")]

print l

pl = [pcr.Primer(k.strip(),v.strip()) for k,v in l]

x = [30.+0.1*i for i in range(0,500)]

for p in pl:
	y = [p.sdss_length(i) for i in x]
	pp.plot(x,y,label=p.name)

pp.legend([p.name for p in pl],loc="upper left")

pp.xlabel("anneal temp")
pp.ylabel("SDSS")
pp.show()