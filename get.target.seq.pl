#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <in.targets.loci.info> <in.seq>\n\n";
	exit 2;
}

my $design = $ARGV[1];

print "#Loci\tSequence\tisOntarget\tRPM\n";
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	##Loci SupportingReads RPM Alignment extra
	my @align = split /;/, $l[3];
	my $seq = '';
	my $ontarget = 'Yes';

	foreach my $i ( 0..$#align ) {
		my $ref = substr( $design, $i, 1);
		my $here = $align[$i];

		if( $ref eq 'N' ) {	## NGG or NNN
			$seq .= $here;
		} else {
			if( $here eq '.' ) {	## perfect match
				$seq .= $ref;
			} elsif( $here eq '-' ) {	## deletion
				$ontarget = 'No';
			} else {
				$here =~ s/^\./$ref/;
				$seq .= $here;
				$ontarget = 'No';
			}
		}
	}
	print join("\t", $l[0], $seq, $ontarget, $l[2]), "\n";
}
close IN;


