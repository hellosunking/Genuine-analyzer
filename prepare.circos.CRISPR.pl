#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;
use Cwd 'abs_path';

if( $#ARGV < 0 ) {
	print STDERR "Usage: $0 <sid> [genome=hg38] [spCas9|spRY]";
	print STDERR "\nPrepare plots for crispr.\n";
	exit 2;
}

my $sid = $ARGV[0];
my $genome = $ARGV[1] || 'hg38';
my $CRISPR_size;	## size of CRISPR sequence
my $enzyme = $ARGV[2] || 'spCas9';
if( $enzyme eq 'spCas9' ) {
	$CRISPR_size = 9288;
} elsif( $enzyme eq 'spRY' ) {
	$CRISPR_size = 8826;
} else {
	print STDERR "ERROR: unknown enzyme!\n";
	exit 1;
}

my $crBinSize   = 5;	## for density plot and links
my $hsBinSize   = 100000;

my $prgDIR = abs_path($0);
$prgDIR =~ s/\/[^\/]+$//;
my $circosDIR = "$prgDIR/circos";

## get genome size and determine the normalization factor for CRISPR
unless( -s "$circosDIR/karyotype.$genome.with.CRISPR.txt" ) {
	print STDERR "ERROR: karyotype file for genome $genome plus CRISPR is missing!\n";
	exit 1;
}

my $genomeSize = 0;
open KAR, "$circosDIR/karyotype.$genome.with.CRISPR.txt";
my $line = <KAR>;
my @l = split /\s/, $line;	##chr - CRISPR 1 0 3088269832 chrC
close KAR;
if( $l[2] eq 'hs0' || $l[2] eq 'mm0') {
	$genomeSize = $l[5];
} else {
	print STDERR "ERROR: karyotype file for genome $genome plus CRISPR is incorrect!\n";
	exit 1;
}
my $CRISPR_norm_factor = int($genomeSize/$CRISPR_size);

## extract target.info
my %targetLoci;
my $onTargetSite = '';
my $normFactor = 0;
my $circleLoop = 0;
open IN, "$sid.targets.loci.info" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	chomp;
	@l = split /\t/;	##chr15:65345197	12016	937.214432727877	.;.;A;.;.;.;A;.;.;.;.;.;.;.;.;.;.;.;.;.;T;G;G	Y	Y
	$targetLoci{$l[0]} = 1;
}
close IN;

## check translocations
my $translocation = '';
my ($hschr, $hspos, $crchr, $crpos);
my (%hsdensity, %crdensity, %links);
open IN, "$sid.translocation" or die( "$!" );
my $metFlag = 0;
while( <IN> ) {
	next if /^#/;
	chomp;
	@l = split /\t/; ##chr19:43730844 chrC:3157 1500
	if( $l[0] =~ /^chrC/ ) {
		if( exists $targetLoci{$l[1]} ) {
			($hschr, $hspos) = split /:/, $l[1];
			($crchr, $crpos) = split /:/, $l[0];
		} else {
			next;
		}
	} elsif( $l[1] =~ /^chrC/ ) {
		if( exists $targetLoci{$l[0]} ) {
			($hschr, $hspos) = split /:/, $l[0];
			($crchr, $crpos) = split /:/, $l[1];
		} else {
			next;
		}
	} else {
		next;
	}

	## density
	my $hsIndex = int($hspos/$hsBinSize)*$hsBinSize;
	$hsdensity{"$hschr:$hsIndex"} += $l[2];
	my $crIndex = int($crpos)/$crBinSize*$crBinSize;
	$crdensity{$crIndex} += $l[2];

	## links
	$links{"$hschr:$hspos-$crIndex"} += $l[2];

	$metFlag = 1;
}
close IN;

unless( $metFlag ) {
	print "WARNING: No CRISPR insertion detected.\n";
	`touch NO_CRISPR_INSERTION`;
	exit;
}

if( $genome =~ /^hg/ ) {
	$crchr = 'hs0';
} elsif( $genome =~ /^mm/ ) {
	$crchr = 'mm0';
}
open OUT, ">$sid.crispr.links" or die( "$!" );
foreach my $loci (sort keys %links) {
	my ($chr, $pos, $crpos) = split /[:-]/, $loci;
	if( $genome =~ /^hg/ ) {
		$chr =~ s/^chr/hs/;
	} elsif( $genome =~ /^mm/ ) {
		$chr =~ s/^chr/mm/;
	}
	my $hsIndex = int($pos/$hsBinSize)*$hsBinSize;
	my $lwd = int(log($links{$loci})) + 1;
	print OUT join("\t", $chr, $hsIndex, $hsIndex+$hsBinSize,
					$crchr, $crpos*$CRISPR_norm_factor, ($crpos+$crBinSize)*$CRISPR_norm_factor,
					"value=$lwd"), "\n";
}
close OUT;

## it seems that density plots are not useful?
open OUT, ">$sid.crispr.density" or die( "$!" );
foreach my $loci (sort keys %hsdensity) {
	my ($chr, $pos) = split /:/, $loci;
	if( $genome =~ /^h/ ) {
		$chr =~ s/^chr/hs/;
	} elsif( $genome =~ /^m/ ) {
		$chr =~ s/^chr/mm/;
	}
	print OUT join("\t", $chr, $pos, $pos+$hsBinSize, log($hsdensity{$loci})), "\n";
}
foreach my $crpos (sort keys %crdensity) {
	print OUT join("\t", $crchr, $crpos*$CRISPR_norm_factor, ($crpos+$crBinSize)*$CRISPR_norm_factor, log($crdensity{$crpos})), "\n";
}
close OUT;

open OUT, ">$sid.circos.crispr.conf" or die( "$!" );
open IN, "$circosDIR/example.CRISPR.conf" or die( "$!" );
while(<IN>) {
	s/PRGPATH/$circosDIR/;
	s/GENOME/$genome/;
	s/SAMPLEID/$sid/;
	print OUT $_;
}
close IN;
close OUT;

