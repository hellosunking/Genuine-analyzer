#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
# v2 add RPM values
# v3 add CRISPR inserts

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <design.sequence> <in.target.loci.info.v2>\n\n";
	exit 2;
}

my $sid = $ARGV[1];
$sid =~ s/.targets.loci.info//;

## HTML header with CSS
print "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"
\"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">
<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\">
<head>
<meta charset=\"UTF-8\">
<title>GENUINE-seq result - $sid</title>
<style type=\"text/css\">
table {
width: 1500px;
table-layout: fixed;
border-collapse: separate;
border-spacing: 0px 2px;
}
td {
font-family: monospace;
font-size: 20px;
text-align: center;
width: 20px;
}
td:nth-child(1) {
text-align: left;
padding-left: 10px;
padding-right: 10px;
width: 160px;
}
td.wider {
padding-left: 5px;
padding-right: 5px;
width: 80px;
}
td.match{color: #606060;}
td.deletion{background-color: #d0d0d0;}

.sA{background-color: #F37777;}
.sC{background-color: #7074B6;}
.sG{background-color: #FAE80C;}
.sT{background-color: #6ABD45;}
.sN{background-color: #d0d0d0;}

tr.withIndel{background-color: grey;}
span.blank {font-size: 12px;}
span.insertion {font-size: 12px;}
</style>
</head>
<body>
<table class=\"GENUINE\">
<tr><td>Loci</td>";

## design sequence
my $design = uc $ARGV[0];
my $len = length $design;

for(my $i=0; $i<$len; ++$i ) {
	print td_letter( substr( $design, $i, 1 ) );
}
print "<td class=\"wider\">Read No.</td><td class=\"wider\">RPM</td><td class=\"wider\">Local circles</td>",
		"<td class=\"wider\">Structural Variants</td><td class=\"wider\">CRISPR insertion</td></tr>\n";
my $ncol = $len + 6;
## TODO: is it necessary to show the "Local" or "Circle" part?

## on/off target sites
open IN, "$ARGV[1]" or die( "$!" );
<IN>;	## header
<IN>;	## header line 2
while( <IN> ) {
	if( /^#/ ) {	## the following loci contains indels, add a separation line here
		print "<tr class=\"withIndel\"><td colspan=$ncol></td></tr>\n";
		next;
	}
	chomp;
	my @l = split /\t/;	###Loci	SupportingReads	RPM Alignment	Circle	Translocation CRISPR extra?
	my @align = split /;/, $l[3];
	my $colored = '';
	foreach my $i ( @align ) {
		my $cell = '';
		if( length $i == 1 ) {	## match or deletion
			$cell = td_letter($i);
		} else {	## insertion
			my $left  = substr($i, 0, 1);
			my $right = substr($i, 1);
			my $blank = '&nbsp;' x length($right);
			$cell = "<td><span class=\"insertion\">$blank</span>";
			$cell .= span_letter($left);
			$cell .= "<span class=\"insertion\">";
			for(my $j=0; $j<length($right); ++$j) {
				$cell .= span_letter( substr($right, $j, 1) );
			}
			$cell .= "</span></td>";
		}
		$colored .= $cell;
	}
	my $rpm = sprintf "%.2f", $l[2];
	print "<tr><td>$l[0]</td>$colored<td>$l[1]</td><td>$rpm</td><td>$l[4]</td><td>$l[5]</td><td>$l[6]</td></tr>\n";
}
close IN;
print "</table><br />\n<hr />Note: Chrome is recommended to view this file.</body>\n</html>";

sub td_letter {
	my $i = shift;

	if( $i eq '.' ) {
		return '<td class="match">&middot;</td>';
	} elsif( $i eq '-' ) {
		return '<td class="deletion">-</td>';
	} else {
		return '<td class="s' . "$i\">$i" . '</td>';
	}
}

sub span_letter {
	my $i = shift;

	if( $i eq '.' ) {
		return '<span class="match">&middot;</span>';
	} elsif( $i eq '-' ) {
		return '<span class="deletion">-</span>';
	} else {
		return '<span class="s' . "$i\">$i" . '</span>';
	}
}
