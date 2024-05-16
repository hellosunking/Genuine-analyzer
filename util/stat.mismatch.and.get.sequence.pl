#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <sid> <target.Seq> [mode=read|loci]\n\n";
	exit 2;
}

my $sid    = $ARGV[0];
my $design = $ARGV[1];
my $mode   = $ARGV[2] || 'read'; 

my (%raw, %norm);
my %sequence;
open IN, "$sid.targets.loci.info" or die( "$!" );
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
	if( $mode eq 'read' ) {
		$raw{$mm}  += $l[1];
	} else {
		$raw{$mm} ++;
	}
	$norm{$mm} += $l[2];
	$sequence{$s} += $l[1];
}
close IN;

#print "Mismatch\t0\t1\t2\t3\t4\t5\n";
my @cnt = map {$raw{$_}||0} (0..5);
print STDERR join("\t", $sid, @cnt ), "\n";

## sequence
foreach my $s ( sort {$sequence{$b} <=> $sequence{$a}} keys %sequence ) {
	print "$s\t$sequence{$s}\n";
}
