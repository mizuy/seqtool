#!/usr/bin/perl -w
# http://www.genenames.org/cgi-bin/hgnc_downloads.cgi
use strict;
use LWP::Simple;
my $url = "http://www.genenames.org/cgi-bin/hgnc_downloads.cgi?".
    "title=HGNC%20output%20data&col=gd_hgnc_id&col=gd_app_sym&col=md_eg_id&col=md_refseq_id&".
    "status=Approved&status=Entry%20Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&".
    "format=text&limit=&submit=submit&.cgifields=&.cgifields=status&".
    ".cgifields=chr&.cgifields=hgnc_dbtag";
my $page = get($url);
print $page;
