#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  : 20230905

use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";
use GENUINEv2 qw/SmithWaterman visualizeSW revcmp loadGenome/;

## v4 deals with Crispr sequence
if( $#ARGV < 3 ) {
	print STDERR "\nUsage: $0 <bg.noise> <design.sequence> <sid> <genome> [max.mismatch=5]\n\n";
	print STDERR "20230907: use Smith-Waterman; RPM and Fold-change to de-noise; allows hotspots in translocation.\n\n";
	exit 2;
}

my $genome = $ARGV[3];

## key parameters
my $minBG  = 2;		## min coverage to be considered as bg noise if provided
my $minFC  = 5;		## min fold-change to noise to be considered
my $minRPM = 0.1;		## TODO: need to tune this parameter
#my $minDirect = 2;	## min for translocation
my $minSupportive = 2;
my $maxSeqDiff = $ARGV[4] || 5;
my $ext  = 60;	## extension of cut point to look for the designed sequence
my $dist = 20;	## distance for checking bg noise

my (%noise, %bgHotSpot);
my $noise_total = 0;
my $lociCNT = 0;
print STDERR "Loading bgNoise: $ARGV[0]\n";
if( $ARGV[0] =~ /\.gz$/ ) {
	open BG, "zcat $ARGV[0] |" or die( "$!" );
} else {
	open BG, "$ARGV[0]" or die( "$!" );
}
while( <BG> ) {
	chomp;
	my @l = split /\t/;	## chrX:YYY count hospot
	my ($c, $p) = split /:/, $l[0];	## chrX:YYY
	if( $#l > 0 ) {	## there is a count
		next if $l[1] < $minBG;
		$noise_total += $l[1];
		$noise{$c}->{$p} = $l[1];
		$bgHotSpot{"$c:$p"} = 1 if $l[2];
	} else {
		$noise{$c}->{$p} = 1;
	}
	++ $lociCNT;
}
close BG;
print STDERR "=> Valid loci: $lociCNT, supporting reads: $noise_total.\n";
$noise_total /= 1e6;	## normalize to 1M reads

my $design = uc $ARGV[1];
my $minScore = length($design) - $maxSeqDiff;
my $g = loadGenome( $genome );

my $sid = $ARGV[2];
my (%loci, %vis, %rawSequence);
my (@withIndel, @noIndel);
## L and R is to avoid the match at head/tail
my $L = '<' x 10;
my $R = '>' x 10;
## cut points: this file CONTAINS the circles and translocations
print STDERR "Loading $sid.cutpoints\n";
my $total_cutpoints = 0;
my %cutpoints;
$lociCNT = 0;
open C, "$sid.cutpoints" or die( "$!" );
while( my $r = <C> ) {
	chomp($r);
	my @l = split /\t/, $r;	##chr6:127625729	34
	next if $l[0]=~/chr[MC]/ || $l[1]<$minSupportive;
	$total_cutpoints += $l[1];
	$cutpoints{$l[0]} = $l[1];
	++ $lociCNT;
}
close C;
print STDERR "=> Valid loci: $lociCNT, supporting reads: $total_cutpoints.\n";
$total_cutpoints /= 1e6;	##normalize to 1M reads

print STDERR "Filtering cutpoints against $design\n";
foreach my $cutlocus ( keys %cutpoints ) {
	## check noise, or this cut point has been reported before
#print STDERR "Dealing $cutlocus\n";
	my ($chr, $pos) = split /:/, $cutlocus;
	my $flag = 0;
	my $rpm_treat = $cutpoints{$cutlocus} / $total_cutpoints;
	next if $rpm_treat < $minRPM;
	foreach my $k ( ($pos-$dist) .. ($pos+$dist) ) {
		if( exists $noise{$chr}->{$k} ) {
			## keep if RPM in treatment is >minFC times higher than background
			if( $noise_total != 0 ) {	## check foldchange
				my $rpm_noise = $noise{$chr}->{$k} / $noise_total;
				if( $rpm_treat / $rpm_noise < $minFC ) {
					$flag = 1;
					last;
				}
			} elsif( exists $loci{"$chr:$k"} ) {	## there is a nearby loci?
				## report a warning?
			} else {	## there is NO counts for noise, so use it qualitatively
				$flag = 1;
				last;
			}
		}
	}
	next if $flag;

	## TODO: deal with Crispr sequence, accept all?
	next if $chr eq 'chrC';

	## check sequence similarity
	my $seq = $L . substr( $g->{$chr}, $pos-$ext, $ext*2 ) . $R;	## check sequence; use M/N's to avoid border-matches
	my $realL = substr( $g->{$chr}, $pos-$ext-10, 10);
	my $realR = substr( $g->{$chr}, $pos+$ext, 10);
	my ($mm1, $os1, $align11, $align12, $ref1, $qry1) = SmithWaterman( $seq, $design );
	my ($mm2, $os2, $align21, $align22, $ref2, $qry2) = SmithWaterman( revcmp($seq), $design );
#	print STDERR "SW: $mm1, $mm2\t$os1, $os2\n";
	my ($aligninfo, $match, $mismatch, $insertion, $deletion);
	my $bestMatch = '';
	if( $mm1 <= $mm2 ) {
		if( $mm1 <= $maxSeqDiff ) {
			$ref1 =~ s/^$L/$realL/;
			$ref1 =~ s/$R$/$realR/;
			($aligninfo, $match, $mismatch, $insertion, $deletion) = visualizeSW( $ref1, $qry1 );
			$bestMatch = $qry1;
		} else {
			next;
		}
	} elsif( $mm2 <= $maxSeqDiff ) {
		## L and R must also be rev-comp
		my $Lrev = revcmp( $realL );
		my $Rrev = revcmp( $realR );
		$ref2 =~ s/$L$/$Lrev/;
		$ref2 =~ s/^$R/$Rrev/;
		($aligninfo, $match, $mismatch, $insertion, $deletion) = visualizeSW( $ref2, $qry2 );
		$bestMatch = $qry2;
	} else {
		next;
	}
	## further filterings
	next if $insertion + $deletion > 1;

	$loci{$cutlocus} = $rpm_treat;
	$vis{$cutlocus}  = $aligninfo;
	$bestMatch=~s/^\s+//;
	$bestMatch=~s/\s+$//;
	$rawSequence{$cutlocus} = $bestMatch;
#	print STDERR "Add loci $l[0]\n";
	if( $insertion || $deletion ) {
		push @withIndel, $cutlocus;
		#print STDERR "ADD $l[0] with Indel ($insertion, $deletion)\n";
	} else {
		push @noIndel, $cutlocus;
		#print STDERR "ADD $l[0] no Indel\n";
	}
}
close C;

my %circleInfo;
print STDERR "Loading $sid.circles\n";
open C, "$sid.circles" or die( "$!" );
#my $header = <C>;
while( my $r = <C> ) {
	next if $r=~/^#/;
	chomp($r);
	my @l = split /\t/, $r;	##chr6:127625729	3
	next if $l[-1] < $minSupportive;
#	print STDERR "Checking $l[0]\n";
	my ($chr, $pos) = split /:/, $l[0];
#	my $cnt = 0;
	foreach my $i ( ($pos-$ext) .. ($pos+$ext) ) {
#		print STDERR "$chr:$i ";
		$circleInfo{"$chr:$i"} = 'Y' if exists $loci{"$chr:$i"};
#		print STDERR "Circle: $l[0] => $chr:$i\n";
	}
#	print STDERR "\n";
}
close C;

my %crispr;
my (%transInfo, %transPair, %transPair_bgHotspot);
print STDERR "Loading $sid.translocation\n";
open T, "$sid.translocation" or die( "$!" );
#$header = <T>;
while( my $r = <T> ) {
	next if $r=~/^#/;
	chomp($r);
	next unless $r=~/\S+/;	## in case of empty line
	my @l = split /\t/, $r; ##chr6:169493400	chrX:107865491	38
	next if $l[-1] < $minSupportive;
#	print STDERR "Checking $l[0], $l[1]\n";

	## check whether both compartments pass
	my $hit_bgHotspot = '';	## hotspots in background noise are ALSO considered
	my ($chr,  $pos ) = split /:/, $l[0];
	my ($chr2, $pos2) = split /:/, $l[1];
	next if $chr eq 'chrC' && $chr2 eq 'chrC';	## both are crispr sequence

	my (@left, @right);
	## TODO: deal with Crispr sequence
	if( $chr ne 'chrC' ) {
		foreach my $i ( ($pos-$ext) .. ($pos+$ext) ) {
			if( exists $loci{"$chr:$i"} ) {
				push @left, "$chr:$i";
			} elsif( exists $bgHotSpot{"$chr:$i"} ) {
				push @left, "$chr:$i";
				$hit_bgHotspot .= 'L';
			}
		}
		next if $#left < 0;
	}

	if( $chr2 ne 'chrC' ) {
		foreach my $i ( ($pos2-$ext) .. ($pos2+$ext) ) {
			if( exists $loci{"$chr2:$i"} ) {
				push @right, "$chr2:$i";
			} elsif ( exists $bgHotSpot{"$chr2:$i"} ) {
				push @right, "$chr2:$i";
				$hit_bgHotspot .= 'R';
			}
		}
		next if $#right < 0;
	}

	if( $chr eq 'chrC' || $chr2 eq 'chrC' ) {
		unless( $hit_bgHotspot=~/L/ ) {
			foreach my $i ( @left ) {
				$crispr{$i} = 'Y';
			}
		}
		unless( $hit_bgHotspot=~/R/ ) {
			foreach my $i ( @right ) {
				$crispr{$i} = 'Y';
			}
		}
	} else {
		next if $hit_bgHotspot=~/L/ && $hit_bgHotspot=~/R/;	## both are bgNoise!
		foreach my $i ( @left, @right ) {
			$transInfo{$i} = 'Y';
		}
		my $k1 = join(",", @left);
		my $k2 = join(",", @right);
		## TODO: check whether k1 and k2 have overlaps!!!
		$transPair{"$k1\t$k2"} += $l[2];
		$transPair_bgHotspot{"$k1\t$k2"} = $hit_bgHotspot || 'N';
		#print STDERR "Adding $k1, $k2\n";
	}
}
close T;

## output the results
open OUT, ">$sid.targets.loci.info" or die("$!");
## SSS: sequence similarity score
#print OUT "#Loci\tSupportingReads\tRPM\tAlignment\tLocal.Circle\tStructural.Variants\tCrispr\tSequence\n";
print OUT "#Loci\tSupportingReads\tRPM\tAlignment\tLocal.Circle\tStructural.Variants\tCrispr\n";
## loci without indels first
print OUT "## Loci without indels\n";
foreach my $i ( sort {$loci{$b} <=> $loci{$a} } @noIndel ) {
	my $r = $loci{$i};
	my $v = $vis{$i};
	my $cnt = $cutpoints{$i};
	print OUT join("\t", $i, $cnt, $r, $v,
		$circleInfo{$i} || 'N', $transInfo{$i} || 'N', $crispr{$i} || 'N'), "\n";
#		$rawSequence{$i}), "\n";
}

if( $#withIndel >= 0 ) {
	print OUT "## Loci with indels\n";
	foreach my $i ( sort {$loci{$b} <=> $loci{$a} } @withIndel ) {
		my $r = $loci{$i};
		my $v = $vis{$i};
		my $cnt = $cutpoints{$i};
		print OUT join("\t", $i, $cnt, $r, $v,
			$circleInfo{$i} || 'N', $transInfo{$i} || 'N', $crispr{$i} || 'N'), "\n";
#			$rawSequence{$i}), "\n";
	}
}
close OUT;

my $tranlocationINFO = '';
foreach my $i ( sort { $transPair{$b} <=> $transPair{$a} } keys %transPair ) {
	my $rpm = $transPair{$i} / $total_cutpoints;
	$tranlocationINFO .= join("\t", $i, $transPair{$i}, $rpm, $transPair_bgHotspot{$i}) . "\n";
}
if( $tranlocationINFO ) {
	open OUT, ">$sid.structural.variants" or die("$!");
	print OUT "#Loci1\tLoci2\tSupportingReads\tRPM\tHitBgHotSpot\n$tranlocationINFO";
	close OUT;
} else {
	`rm -f $sid.structural.variants && touch $sid.structural.variants`;
}

sub calc_diff {
	my $ref = uc shift;
	my $qry = uc shift;

	my $rsize = length( $ref );
	my $qsize = length( $qry );
	my $minMatch = $qsize >> 1;

	if( $rsize < $qsize ) {
		print STDERR "ERROR: Query is longer than reference!\n";
		return 100;
	}

	my ($bestScore, $bestIndex) = (-100, -100);
	foreach my $i ( 0 .. $rsize ) {
		my $check = substr( $ref, $i, $qsize );
		last if length($check) < $minMatch;

		my $score = calc_align_score( $check, $qry );
		if( $score > $bestScore ) {
			$bestScore = $score;
			$bestIndex = $i;
		}
	}

	return $bestScore;
}

sub calc_align_score {
	my $a = shift;
	my $b = shift;
	my $s = (length($a) < length($b)) ? length($a) : length($b);

	my $v = 0;
	for(my $j=0; $j!=$s; ++$j) {
		my $aa = substr($a, $j, 1);
		my $bb = substr($b, $j, 1);
		++ $v if $aa eq $bb || $aa eq 'N' || $bb eq 'N';
	}

	return $v;
}

sub isOntarget {
	my $qry = shift;
	my $ref = shift;

	return "N" if length($qry) != length($ref);
	my $len = length($ref);
	for( my $i=0; $i<$len; ++$i ) {
		if( substr($ref, $i, 1) != 'N' ) {
			return 'N' if substr($ref, $i, 1) ne substr($qry, $i, 1);
		}
	}
	return 'Y';
}

