all: ucsc.tab hgnc.tab

ucsc.tab:
	curl --output ucsc.tab "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=310886027&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=refGene&hgta_table=refGene&hgta_regionType=genome&position=chr21%3A33031597-33041570&hgta_outputType=primaryTable&hgta_outFileName=&submit=submit&hgta_doTopSubmit=1"



hgnc.tab:
	curl --output hgnc.tab "http://www.genenames.org/cgi-bin/hgnc_downloads?col=gd_app_sym&col=md_eg_id&status=Approved&status=Entry%20Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit"
