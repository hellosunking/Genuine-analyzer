#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <targets.loci.info> <target.Seq>\n\n";
	exit 2;
}

my $design = $ARGV[1];

my (%raw, %norm);
my %sequence;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	##chr1:185087641	773	93.5960561309957	.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;A;.;.	Y	N

	my $mm = 0;
	my @align = split /;/, $l[3];
	my $s = '';
	foreach my $i ( 0..$#align ) {
		my $ref = substr( $design, $i, 1 );
		if( $align[$i] eq '.' ) {
			$s .= $ref;
		} else {
			++ $mm if $ref ne 'N';
			if( $align[$i] eq '-' ) {	## deletion
			} elsif( $align[$i] =~ s/^\.// ) {	## insertion
				$s .= $ref . $align[$i];
			} elsif( $align[$i] =~ /^[ACGT]$/i ) {
				$s .= $align[$i];
			} else {
				print STDERR "ERROR $l[3]\n";
			}
		}
	}

	next unless $mm;	## ignore on-target site
	my ($c,$p) = split /:/, $l[0];
	print join("\t", $c, $p-1, $p) . "\n";

#	my $bed = join("\t", $c, $p-100, $p+100) . "\n";
#	print $bed x $l[1];
}
close IN;

