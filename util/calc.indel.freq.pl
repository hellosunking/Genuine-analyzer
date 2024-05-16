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
my $sid = shift;

my %targets;
open IN, "$sid.targets.loci.info" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;
	$targets{$l[0]} = 1;
}
close IN;

my %dist;
open IN, "zcat unc.pass.cutpoints.gz ext.pass.cutpoints.gz |" or die( "$!" );
while( <IN> ) {
	my @l = split /\t/;	##Pair	chr20	32612679	-	chr20	32612681	+	A00808:1380:HV35HDSX5:4:1126:14633:32503#-1#72
	next unless $l[0] eq 'Pair' && $l[1] eq $l[4];

	## check whether it covers a valid target site
	my $hitTarget = 0;
	for( my $i=$l[2]-$ext; $i<=$l[2]+$ext; ++$i ) {
		if( exists $targets{"$l[1]:$i"} ) {
			$hitTarget = 1;
			last;
		}
	}
	next unless $hitTarget;
	$hitTarget = 0;
	for( my $i=$l[5]-$ext; $i<=$l[5]+$ext; ++$i ) {
		if( exists $targets{"$l[1]:$i"} ) {
			$hitTarget = 1;
			last;
		}
	}
	next unless $hitTarget; 

	## check distance
	my $d = $l[5] - $l[2];
	if( $d == 0 ) {	## 1-bp insertion !
		$dist{0} ++;
	} elsif( $d < 100 && $d > -100 ) {	## small indel, check direction
		if( $l[3] eq '-' && $l[6] eq '+' ) {
			$dist{$d} ++;
		} elsif( $l[3] eq '+' && $l[6] eq '-' ) {
			$dist{-$d} ++;
		} else {
			print STDERR "ERROR DIST: $l[-1]";
		}
	} else {	## large deletion, do not check direction
		$dist{1000} ++;
	}
}
close IN;

print "#$sid\n";
my @i = sort {$a<=>$b} keys %dist;
print join("\t", @i), "\n";
my @j = map {$dist{$_}} @i;
print join("\t", @j), "\n";

