#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <sid> [rmdup=no|yes]\n\n";
	exit 2;
}

my $sid = $ARGV[0];
my $rmdup = $ARGV[1] || 'no';

## preprocessing
my %trim;
open IN, "$sid.trim.log" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##Total   9814563
	$trim{$l[0]} = $l[1];
}
close IN;

open IN, "$sid.flash.log" or die( "$!" );
while( <IN> ) {
	chomp;
	if( /Combined pairs:\s+(\S+)/ ) {
		$trim{'FLASH'} = $1;
		last;
	}
}
close IN;

## spike-in
my ($spikeinERR, $spikeinPass) = ( 0, 0 );
open IN, "cat ext.spikeIn.stat unc.spikeIn.stat |" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/; ##Border	1570746
	if( $l[0] eq 'NoSpikeIn' || $l[0] eq 'Discard' ) {
		$spikeinERR += $l[1];
	} else {
		$spikeinPass += $l[1];
	}
}
close IN;

## De-dup, optional
my $dedup = 0;
if( $rmdup ne 'no' ) {
	open IN, "cat *.uniq.stat |" or die( "$!" );
	while( <IN> ) {
		chomp;
		my @l = split /\t/; ##Uniq    243786
		$dedup += $l[1] if $l[0] eq 'Uniq';
	}
	close IN;
} else {
	$dedup = $spikeinPass;
}

## alignment and cut-point extraction
my %align;
open IN, "cat *.cutpoints.stat|" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/; ##Uniq    243786
	$align{$l[0]} += $l[1];
}
close IN;

## targets and translocations
my ($targets, $trans) = ( 0, 0 );
open IN, "$sid.targets.loci.info" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	++ $targets;
}
close IN;

open IN, "$sid.structural.variants" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	++ $trans;
}
close IN;

## statistics
my $preprocessed = $trim{Total}-$trim{Dropped};

print "Sid\t$sid",
"\nTotal\t", $trim{Total}, "\t100%",
"\nPreprocessed\t", $trim{Total}-$trim{Dropped}, "\t", sprintf("%.2f%%", 100-$trim{Dropped}/$trim{Total}*100),
"\n  Stitch-able\t", $trim{FLASH}, "\t", sprintf("%.2f%%", $trim{FLASH}/$preprocessed*100),
"\nSpike-in\t", $spikeinPass, "\t", sprintf("%.2f%%", $spikeinPass/$preprocessed*100);

print "\nUnique\t$dedup\t", sprintf("%.2f%%", $dedup/$trim{FLASH}*100) if $rmdup ne 'no';
print
"\nAligned\t", $align{Alignable}, "\t", sprintf("%.2f%%", $align{Alignable}/$dedup*100),
"\nFiltered\t", $align{Pass}, "\t", sprintf("%.2f%%", $align{Pass}/$align{Alignable}*100),
"\nCutting loci\t$targets",
"\n  Structual variants\t$trans\n";

