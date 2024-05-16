#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";
use GENUINEv2 qw/revcmp $spikeSeq $spikeSize/;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <R1.fq> <R2.fq> <out.prefix> [max.err=2]\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $minSpikeLen = 10;	## request at least 10 basepairs matching the spike in sequence
my $minSpikeErr = $ARGV[3] || 2;	## allows 1 mismatch
#my $minFragSize = 36;
## dual index
my $index1 = substr( $spikeSeq, 0, 3 );
my $index2 = substr( $spikeSeq, 3, 3 );

## NOTE: if both R1 and R2 contains part of the spike-in sequence in tail, should we keep it?
my $minMappable = 16;	## minmum cycles to be considered as mappable

open BRD, ">$ARGV[2].border.fq" or die( "$!" );
open OR1, ">$ARGV[2].pass.R1.fq" or die( "$!" );
open OR2, ">$ARGV[2].pass.R2.fq" or die( "$!" );

my %stat;
open R1, "$ARGV[0]" or die( "$!" );
open R2, "$ARGV[1]" or die( "$!" );
while( my $rid = <R1> ) {
	chomp( $rid );
	$rid =~ s/\s.*//;	## remove spaces because we will add labels to readID
## Read 1
	my $s1 = <R1>;
	chomp($s1);
	<R1>;
	my $q1 = <R1>;
	chomp($q1);

## Read 2
	<R2>;
	my $s2 = <R2>;
	chomp($s2);
	<R2>;
	my $q2 = <R2>;
	chomp($q2);

	my ($pos1, $slabel1, $fq1) = dealSpikeIn( $s1, $q1 );
	my ($pos2, $slabel2, $fq2) = dealSpikeIn( $s2, $q2 );
#	print STDERR "$rid\t$slabel1 ($pos1)\t$slabel2 ($pos2)\n";

#	++ $stat{"Total"};
	if( $slabel1 eq 'Head' ) {	## keep R2 only
		if( $slabel2 eq 'NoSpikeIn' ) {
			print BRD "$rid#$pos1#$pos2\n$fq2";
			++ $stat{"Border"};
		} else {	## invalid
			++ $stat{"Discard"};
		}
	} elsif( $slabel1 eq 'Pass' ) {
		if( $slabel2 eq 'NoSpikeIn' || $slabel2 eq 'Pass' ) {
			print OR1 "$rid#$pos1#$pos2\n$fq1";
			print OR2 "$rid#$pos1#$pos2\n$fq2";
			++ $stat{"Pass"};
		} else {
			++ $stat{"Discard"};
		}
	} else {	## R1 does not contain spike-in
		if( $slabel2 eq 'Head' ) {	## keep R1 only
			print BRD "$rid#$pos1#$pos2\n$fq1";
			++ $stat{"Border"};
		} elsif ( $slabel2 eq 'Pass' ) {
			print OR1 "$rid#$pos1#$pos2\n$fq1";
			print OR2 "$rid#$pos1#$pos2\n$fq2";
			++ $stat{"Pass"};
		} else {
			## discard, both R1 and R2 does not contain spike-in
			++ $stat{"Discard"};
		}
	}
	## note: I only consider NoSpikeIn + NoSpikeIn/Pass/Head, which covers >95% of the reads
	## others like Pass + Pass is unreasonable and very few
	## TODO: do I need to change NoSpikeIn+head to head+NoSpikeIn similar to SE data?
}
close R1;
close R2;
close BRD;
close OR1;
close OR2;

open STAT, ">$ARGV[2].spikeIn.stat" or die("$!");
foreach my $k ( sort keys %stat ) {
	print STAT "$k\t$stat{$k}\n";
}
close STAT;

sub dealSpikeIn {
	my $s = shift;
	my $q = shift;

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

	my $statLabel = 'NoSpikeIn';
	my $slen = length( $s );
	foreach $i ( @uniqSeed ) {
		last if $i + $spikeSize > $slen;

		my $e = checkSpikeIn( substr($s, $i, $spikeSize), $spikeSeq );
		if( $e != 0 )	{ ## found the spike in, update the sequence
			my ($ss, $qq) = ( '', '' );

			if( $i < $minMappable ) {
				## ss and qq are meaningless as I won't use them in this case
#				$ss = substr( $s, $i+$spikeSize );
#				$qq = substr( $q, $i+$spikeSize );
				$statLabel = 'Head';
			} else {
				$ss = substr( $s, 0, $i );
				$qq = substr( $q, 0, $i );
				## TODO: discard the tailing cycles, which could be short and unmappable
				$statLabel = 'Pass';
			}
			return ($i, $statLabel, "$ss\n+\n$qq\n");
		}
	}

	## spike-in not found
	return (-1, $statLabel, "$s\n+\n$q\n");
}

sub checkSpikeIn {
	my $a = shift;
	my $b = shift;

	my $len = length($a);
	$len = length($b) if length($b) < $len;
	return 0 if $len < $minSpikeLen;

	my $d = 0;
	for(my $i=0; $i!=$len; ++$i) {
		++ $d if substr($a, $i, 1) ne substr($b, $i, 1);
		return 0 if $d > $minSpikeErr;
	}

	return 1;
}
