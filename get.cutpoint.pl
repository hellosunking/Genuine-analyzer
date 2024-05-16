#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <in.cutpoints> <out.prefix>\n\n";
	exit 2;
}

my $circleDist = 100;	## TODO: consider 50?
my $distance   = 100;	## pool nearby cut points
my $minTransDist = 1000;	## min distance to call a putative translocation

my %cutpoints;
my %circle;	## cut points that are perfectly repaired
my %trans;	## potential translocations

## TODO: infer the "designed" loop
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;

	if( $l[0] eq 'Cutpoint' ) {	#Cutpoint	chr15	66442522	+	A00917:1167:HMCC7DSX5:3:1101:5276:1078
		$cutpoints{$l[1]}->{$l[2]} ++;
	} else {	##Pair	chrX	152586325	-	chr1	186197872	+	A00917:1167:HMCC7DSX5:3:1101:29686:5274
		if( $l[1] eq $l[4] && abs($l[2]-$l[5])<$circleDist ) {	## NOTE: consider strand?
			$circle{$l[1]}->{$l[2]} ++;
			$cutpoints{$l[1]}->{$l[2]} ++;
		} else {
			$cutpoints{$l[1]}->{$l[2]} ++;
			$cutpoints{$l[4]}->{$l[5]} ++;

			next if $l[1] eq $l[4] && abs($l[2]-$l[5])<$minTransDist;	## consider two cutpoints, but NOT translocations
		
			## potential translocations, record the one with smaller chr/coordinate as left-side
			if( $l[1] lt $l[4] ) {
				$trans{"$l[1]:$l[2]\t$l[4]:$l[5]"} ++;
			} elsif ( $l[1] gt $l[4] ) {
				$trans{"$l[4]:$l[5]\t$l[1]:$l[2]"} ++;
			} else {
				if( $l[2] < $l[5] ) {
					$trans{"$l[1]:$l[2]\t$l[4]:$l[5]"} ++;
				} else {
					$trans{"$l[4]:$l[5]\t$l[1]:$l[2]"} ++;
				}
			}
		}
	}
}
close IN;

pool_nearby_loci( \%cutpoints, "$ARGV[1].cutpoints" );
pool_nearby_loci( \%circle, "$ARGV[1].circles" );
deduce_trans( \%trans, "$ARGV[1].translocation" );

sub pool_nearby_loci {
	my $cut = shift;
	my $filename = shift;

	open OUT, ">$filename" or die( "$!" );
	foreach my $chr ( keys %$cut ) {
		my $here = $cut->{$chr};
		my @loci = sort { $here->{$b} <=> $here->{$a} } keys %{$here};

		my %cnt;
		$cnt{$loci[0]} = $here->{$loci[0]};
		## for each loci, pool it to its nearby one that has a larger hits
		for( my $i=1; $i<=$#loci; ++$i ) {
			my $j;
			for( $j=0; $j<$i; ++$j ) {
				next if $here->{$loci[$j]} <= 0;	## this site has been pooled
				if( abs($loci[$i] - $loci[$j]) <= $distance ) {	## found, add this one to that loci, and mask this site
					$cnt{$loci[$j]} += $here->{$loci[$i]};
					$here->{$loci[$i]} = -1;
#					print STDERR "$chr:$loci[$i] is pooled to $loci[$j]\n";
					last;
				}
			}
			if( $j >= $i ) {	## not found, consider this site as a new loci
				$cnt{$loci[$i]} = $here->{$loci[$i]};
#				print STDERR "Add $chr:$loci[$i]\n";
			}
		}

		foreach my $k ( keys %cnt ) {
			print OUT "$chr:$k\t", $cnt{$k}, "\n";
		}
	}
	close OUT;
}

sub deduce_trans {
	my $trans = shift;
	my $filename = shift;

	my @loci = sort { $trans->{$b} <=> $trans->{$a} } keys %{$trans};
	open OUT, ">$filename" or die( "$!" );

	if( $#loci < 0 ) {	## no record
		close OUT;
		return;
	}

	my %cnt;
	$cnt{ $loci[0] } = $trans->{$loci[0]};
	for( my $i=1; $i<=$#loci; ++$i ) {
		my @qry = split /[:\t]/, $loci[$i];
		my @cand = sort { $cnt{$b} <=> $cnt{$a} } keys %cnt;

		my $flag = 0;	## whether this loci is pool-able
		foreach my $target ( @cand ) {
			my @ref = split /[:\t]/, $target;
			if( $qry[0] eq $ref[0] && abs($qry[1]-$ref[1])<=$distance && $qry[2] eq $ref[2] && abs($qry[3]-$ref[3])<=$distance ) {	## found
				$cnt{ $target } += $trans->{$loci[$i]};
				$flag = 1;
#				print STDERR "$loci[$i] is pooled to $target\n";
			}
		}
		unless( $flag ) {	## could not be pooled with current loci
			$cnt{ $loci[$i] } = $trans->{$loci[$i]};
		}
	}

	foreach my $i ( keys %cnt ) {
		print OUT "$i\t", $cnt{$i}, "\n";
	}
	close OUT;
}

