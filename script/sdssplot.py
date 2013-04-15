import seqtool.pcr as pcr
import matplotlib.pyplot as pp

l = [x.split(":") for x in """
MA:  CGAACGACGCCCTACGAAAACG
MS:  CGTCGCGGGTCGGAGTTTC
UA:  CAACAAACCAAACAACACCCTACAAA
US:  GGAGGTGTTGTGGGTTGGAGTTTTG
MS2: CGCGGCGGTCGTAGTC
MA3: GACGACCGCGAAATCGACG
US2: GAGTGGGTGTGGTGGTTGTAGTTG
UA3: AAAAACAACAACAACAACAACCACAA
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
