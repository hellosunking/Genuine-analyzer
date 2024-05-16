#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";
use GENUINEv2 qw/revcmp $spikeSeq $spikeSize/;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <in.fq> <out.prefix> [max.err=2]\n\n";
	exit 2;
}

my $spikeErr  = $ARGV[2] || 2;	##  mismatches allowed in SpikeIn
my $minFragSize = 36;

## dual index
my $index1 = substr( $spikeSeq, 0, 3 );
my $index2 = substr( $spikeSeq, 3, 3 );

my $minMappable = 16;	## minmum cycles to be considered as potentially mappable
open BD, ">$ARGV[1].border.fq" or die( "$!" );
open R1, ">$ARGV[1].pass.R1.fq" or die( "$!" );
open R2, ">$ARGV[1].pass.R2.fq" or die( "$!" );

my %stat;
open IN, "$ARGV[0]" or die( "$!" );
while( my $rid = <IN> ) {
	chomp( $rid );
	$rid =~ s/\s.*//;	## remove spaces because we will add labels to readID
	my $s = <IN>;
	chomp($s);
	<IN>;
	my $q = <IN>;
	chomp($q);

	$s = uc $s;
	my %seed;
	my $i = 0;
	while( 1 ) {
		$i = index( $s, $index1, $i );
		last if $i < 0;
		$seed{$i} = 1;
		++ $i;
	}
	$i = 3;
	while( 1 ) {
		$i = index( $s, $index2, $i );
		last if $i < 0;
		my $j = $i - 3;
		$seed{$j} = 1;
		++ $i;
	}
	my @uniqSeed = sort {$a<=>$b} keys %seed;

	my $label = 'NoSpikeIn';
	my $slen = length( $s );
	my $hasSpikeIn = 0;
	my ($s1, $s2, $q1, $q2);
	foreach $i ( @uniqSeed ) {
		last if $i+$spikeSize > $slen+$spikeErr;

		$hasSpikeIn = checkSpikeIn( substr($s, $i, $spikeSize), $spikeSeq );
		if( $hasSpikeIn )	{ ## hit the spike in
#			print STDERR "$rid Spike-in found at $i\n";
			if( $i < $minMappable ) {	## head
				$label = 'Head';
				$stat{'Head'} ++;
				## do rev-comp for heads as we are interested in the end faraway from the spike-in
				$s1 = revcmp( substr($s, $i+$spikeSize) );
				$q1 = revcmp( substr($q, $i+$spikeSize) );

				## the $i $slen are used to aid duplication removal
				## because CRISPR always cut the same site, leading to identical read start/end in Genuine-seq,
				## but the insert size could be different
				print BD "$rid#$i#$slen\n$s1\n+\n$q1\n";
			} elsif( $i+$spikeSize > $slen-$minMappable ) {	## tails
				$label = 'Tail';
				$stat{'Tail'} ++;
				$s1 = substr( $s, 0, $i );
				$q1 = substr( $q, 0, $i );

				print BD "$rid#$i#$slen\n$s1\n+\n$q1\n";
			} else {
				$label = 'Pass';
				$stat{'Pass'} ++;
				$s1 = substr( $s, 0, $i );
				$q1 = substr( $q, 0, $i );
				$s2 = revcmp( substr($s, $i+$spikeSize) );
				$q2 = revcmp( substr($q, $i+$spikeSize) );

				print R1 "$rid#$i#$slen\n$s1\n+\n$q1\n";
				print R2 "$rid#$i#$slen\n$s2\n+\n$q2\n";
			}
			last;
		}
	}
	++ $stat{'NoSpikeIn'} unless $hasSpikeIn;
}
close IN;
close BD;
close R1;
close R2;

open STAT, ">$ARGV[1].spikeIn.stat" or die("$!");
foreach my $k ( sort keys %stat ) {
	print STAT "$k\t$stat{$k}\n";
}
close STAT;

sub checkSpikeIn {
	my $a = shift;
	my $b = shift;

	my $len = length($a);
	$len = length($b) if length($b) < $len;
	return 0 if $len < $spikeSize;

	my $d = 0;
	for(my $j=0; $j!=$len; ++$j) {
		++ $d if substr($a, $j, 1) ne substr($b, $j, 1);
		return 0 if $d > $spikeErr;
	}
	return 1;
}

