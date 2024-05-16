#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <sid>\n\n";
	exit 2;
}

my $ext = 100;
my $minDist = 10000;

my (%loci, %pair);
my $flag = 0;
open IN, "$ARGV[0].structural.variants" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##chr9:96281861	chr9:123375901	1994	150.331407623687	N
	next unless $l[0] =~ /(chr\S+):(\d+)/;
	my ($a, $p1) = ($1, $2);
	next unless $l[1] =~ /(chr\S+):(\d+)/;
	my ($b, $p2) = ($1, $2);
	next unless $a eq $b;
	$loci{$l[0]} = 1;
	$loci{$l[1]} = 1;

	if( $p1 < $p2 ) {
		$pair{"$a:$p1-$p2"} = 1;
	} else {
		$pair{"$a:$p2-$p1"} = 1;
	}
	$flag = 1;
}
close IN;

unless( $flag ) {
	print STDERR "WARNING: No candidate found!\n";
	exit;
}

my %count;
open IN, "zcat ext.pass.cutpoints.gz unc.pass.cutpoints.gz |" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/; ##Pair	chrX	129906664	-	chrX	72308233	+	A00917:1266:H5KMLDSX7
	next unless $l[0] eq 'Pair' && $l[1] eq $l[4] && abs($l[2]-$l[5]) >= $minDist;
	my $chr = $l[1];
	my ($p1, $s1, $p2, $s2);
	if( $l[2] < $l[5] ) {
		($p1, $s1, $p2, $s2) = ($l[2], $l[3], $l[5], $l[6]);
	} else {
		($p1, $s1, $p2, $s2) = ($l[5], $l[6], $l[2], $l[3]);
	}
	my (@left, @right);
	foreach my $i ( ($p1-$ext) .. ($p1+$ext) ) {
		push @left,  $i if exists $loci{"$chr:$i"};
	}
	foreach my $j ( ($p2-$ext) .. ($p2+$ext) ) {
		push @right, $j if exists $loci{"$chr:$j"};
	}
	next unless $#left>=0 && $#right>=0;
	if($#left>=1 || $#right>=1) {
		print STDERR "Too many hits: ", join(",", @left),  " -- ", join(",", @right),  "\n";
	}
	foreach my $L ( @left ) {
		foreach my $R ( @right ) {
			next unless exists $pair{"$chr:$L-$R"};
			++ $count{"$chr:$L-$R"}->{"$s1$s2"};
		}
	}
}
close IN;

my @strands = ('++', '+-', '-+', '--');
print join("\t", "Loci", @strands), "\n";
foreach my $loci (sort keys %count) {
	my $c = $count{$loci};
	my @n = map {$c->{$_} ||0} @strands;
	print join("\t", $loci, @n), "\n";
}

