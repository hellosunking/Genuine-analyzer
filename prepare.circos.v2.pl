#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;
use Cwd 'abs_path';

if( $#ARGV < 0 ) {
	print STDERR "Usage: $0 <sid> <design.seq> [genome=hg38]";
	print STDERR "\nv2 does not plot self-circles.\n";
	exit 2;
}

my $genome = $ARGV[2] || 'hg38';
$genome =~ s/sp.*//;	## remove the spCase/spRY suffix

my $prgDIR = abs_path($0);
$prgDIR =~ s/[^\/]+$//;
my $circosDIR = "$prgDIR/circos";

my $circleDist = 100;   ## TODO: consider 50?
my $bin = 1e5;
my $enlarge_factor = 3;
## color scheme
my $color_ontarget_circle  = "vdgreen";
my $color_offtarget_circle = "lgreen";

my $color_on_off_trans      = "blue";
my $color_on_hotspot_trans  = "dblue";
my $color_off_off_trans     = "red";
my $color_off_hotspot_trans = "dgrey";

my $sid = $ARGV[0];
my $design = $ARGV[1];

## extract target.info
my %targetLoci;
my $onTargetSite = '';
my $normFactor = 0;
my $circleLoop = 0;
open IN, "$sid.targets.loci.info" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	##chr15:65345197	12016	937.214432727877	.;.;A;.;.;.;A;.;.;.;.;.;.;.;.;.;.;.;.;.;T;G;G	Y	Y
	$normFactor = $l[1]/$l[2] if $normFactor<0.001;

	## find the ontarget site
	my @align = split /;/, $l[3];
	my $mismatch = 0;
	foreach my $i ( 0..$#align ) {
		if( $align[$i] ne '.' ) {
			if( substr( $design, $i, 1 ) ne 'N' ) {
				$mismatch = 1;
				last;
			}
		}
	}

	my ($c, $p) = split /:/, $l[0];
	if( $genome eq 'hg38' ) {
		$c =~ s/chr/hs/;
	} elsif( $genome eq 'mm10' ) {
		$c =~ s/chr/mm/;
	}
	$p = int($p/$bin) * $bin;
	my $k = join("\t", $c, $p, $p+$bin);

	if( $mismatch == 0 ) {	## ontarget site
		$circleLoop = 1;
		$onTargetSite = $l[0];

		## always plot the on-target site even it does not have any circles/translocations
		if( $l[4] eq 'Y' ) {
			$targetLoci{$k} = "OnTargetCircle";
		} else {
			$targetLoci{$k} = "OnTarget";
		}
	} elsif( $l[4] eq 'Y' ) {	## check circles
		$circleLoop = 1;
		$targetLoci{$k} = "OffTargetCircle";
	}
}
close IN;

#if( $onTargetSite ) {
#	print STDERR "INFOR: on-target site is $onTargetSite.\n";
#} else {
#	print STDERR "WARNING: No on-target site!\n";
#}
#print STDERR "INFOR: Normalization factor is $normFactor.\n";

## check translocations
my $translocation = '';
open IN, "$sid.structural.variants" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/; ##chr15:67840045	chr9:130998106	1058	82.5210444262728	L
	next if $l[4]=~/L/ && $l[4]=~/R/;	## both are hotspots

	## check whether left/right contains ontarget sites
	my $isOntarget = 0;
	foreach my $loci ( split /,/, $l[0] ) {
		if( $loci eq $onTargetSite ) {
			$isOntarget = 1;
			last;
		}
	}
	foreach my $loci ( split /,/, $l[1] ) {
		if( $loci eq $onTargetSite ) {
			$isOntarget = 2;
			last;
		}
	}

	## get loci
	$l[0] =~ s/,.*//;
	my ($c1, $p1) = split /:/, $l[0];
	if( $genome eq 'hg38' ) {
		$c1 =~ s/chr/hs/;
	} elsif( $genome eq 'mm10' ) {
		$c1 =~ s/chr/mm/;
	}
	$p1 = int($p1/$bin) * $bin;
	$l[1] =~ s/,.*//;
	my ($c2, $p2) = split /:/, $l[1];
	if( $genome eq 'hg38' ) {
		$c2 =~ s/chr/hs/;
	} elsif( $genome eq 'mm10' ) {
		$c2 =~ s/chr/mm/;
	}
	$p2 = int($p2/$bin) * $bin;

	## get translocation type for color
	my $color_here;
	if( $l[4]=~/[LR]/ ) {
		if( $isOntarget ) {
			$color_here = $color_on_hotspot_trans;
		} else {
			$color_here = $color_off_hotspot_trans;
		}
	} else {
		if( $isOntarget ) {
			$color_here = $color_on_off_trans;
		} else {
			$color_here = $color_off_off_trans;
		}
	}

	my $t = int(log($l[3])) * $enlarge_factor;
	$t = 1 if $t < 1;
	$translocation .= join("\t", $c1, $p1, $p1+$bin, $c2, $p2, $p2+$bin, "color=$color_here,value=$t") . "\n";
#	print STDERR "Translocation: $l[0] - $l[1]\n";

	my ($left_type, $right_type) = ('', '');
	if( $l[4]=~/L/ ) {
		$left_type  = "HotSpot";
		$right_type = ($isOntarget) ? "OnTarget":"OffTarget";
	} elsif( $l[4]=~/R/ ) {
		$right_type  = "HotSpot";
		$left_type = ($isOntarget) ? "OnTarget":"OffTarget";
	} else {
		if( $isOntarget == 0 ) {
			$left_type  = "OffTarget";
			$right_type = "OffTarget";
		} elsif( $isOntarget == 1 ) {
			$left_type  = "OnTarget";
			$right_type = "OffTarget";
		} else {
			$left_type  = "OffTarget";
			$right_type = "OnTarget";
		}
	}
	my $k = join("\t", $c1, $p1, $p1+$bin);
	$targetLoci{$k} = $left_type  unless exists $targetLoci{$k};
	$k = join("\t", $c2, $p2, $p2+$bin);
	$targetLoci{$k} = $right_type unless exists $targetLoci{$k};
}
close IN;

if( $circleLoop || $translocation ) {
	open OUT, ">$sid.links" or die( "$!" );
	print OUT $translocation;
	close OUT;

	open OUT, ">$sid.loci" or die( "$!" );
	foreach my $i ( sort keys %targetLoci ) {
		print OUT "$i\t0\ttype=", $targetLoci{$i}, "\n";
	}
	close OUT;

	open OUT, ">$sid.circos.conf" or die( "$!" );
	open IN, "$circosDIR/exampleV2.circos.conf" or die( "$!" );
	while(<IN>) {
		s/PRGPATH/$circosDIR/;
		s/GENOME/$genome/;
		s/SAMPLEID/$sid/;
		print OUT $_;
	}
	close IN;
	close OUT;
} else {
	print STDERR "WARNING: No loops or translocations found! Circos plot will be ommited.\n";
}

