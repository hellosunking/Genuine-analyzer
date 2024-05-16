#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <structural.variants> <ext.pass.cutpoints> <unc.pass.cutpoints>\n\n";
	exit 2;
}

my $ext = 50;    ## extensions from cut points

my %hits;
my %pair;

open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	##chr1:98882095	chr12:1878905	63
	$hits{$l[0]} = 1;
	$hits{$l[1]} = 1;

	$pair{"$l[0]\t$l[1]"} = 1;
}
close IN;

open IN, "less $ARGV[1] $ARGV[2] |" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/; ##Pair	chr14	93073272	-	chr3	80758955	+	A00917:1167:HMCC7DSX5:3:1101:1642:1110#97#235
	next unless $l[0] eq 'Pair';

	my (%left, %right);
	my ($flag_left, $flag_right) = ( 0, 0 );
	for(my $i=$l[2]-$ext; $i<=$l[2]+$ext; ++$i ) {
		if( exists $hits{"$l[1]:$i"} ) {
			$left{"$l[1]:$i"} = 1;
			$flag_left = 1;
		}
	}
	next unless $flag_left;
	for(my $i=$l[5]-$ext; $i<=$l[5]+$ext; ++$i ) {
		if( exists $hits{"$l[4]:$i"} ) {
			$right{"$l[4]:$i"} = 1;
			$flag_right = 1;
		}
	}
	next unless $flag_right;

	foreach my $i ( keys %left ) {
		foreach my $j ( keys %right ) {
			if( exists $pair{"$i\t$j"} ) {
				print "$i\t$j\t$l[-1]\n";
			} elsif( exists $pair{"$j\t$i"} ) {
				print "$j\t$i\t$l[-1]\n";
			}
		}
	}
}
close IN;

