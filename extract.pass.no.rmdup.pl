#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

use FindBin qw/$Bin/;
use lib "$Bin";
use GENUINEv2 qw/get_clip_size_from_cigar extract_cutpoint $min_mapQ $min_mapped_ratio $max_unmapped_allowed/;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <pass.bam> <out.prefix>\n\n";
	exit 2;
}

my $rescueMode = 1;
my $distance   = 100;   ## if the cut points are close to each other, consider it as a circle

my %met;
my ($alignable, $discarded, $pass) = (0,0,0);
my $result = '';

## start analysis
my $lastID = '';
my @lastRead = ();
if( $ARGV[0] =~ /bam$/ ) {
	open IN, "samtools view $ARGV[0] |" or die( "$!" );
} else {
	open IN, "$ARGV[0]" or die( "$!" );
}
my $cnt = 0;
while( <IN> ) {
	next if /^@/;
	chomp;
	my @l = split /\t/;	##A00133:521:H3TMGDSX3:1:1104:8467:14372	0	chr5	149462986	255	140M122S	*	0
	next if $l[4] < $min_mapQ;	## low quality alignment
	next if ($l[1] & 256) || ($l[1] & 512) || ($l[1] & 1024);	## 2nd alignment, or QC failed, or PCR duplicate
	my @key = ($l[1], $l[2], $l[3], $l[4], $l[5]);

	if( $l[0] ne $lastID ) {	## meet a new read, analyze the previous one
		analyze_segment( $lastID, \@lastRead );
		$lastID = $l[0];
		@lastRead = ();
		push @lastRead, \@key;
		++ $cnt;
	} else {	## same read to the previous record
		push @lastRead, \@key;
	}
}
close IN;
analyze_segment( $lastID, \@lastRead );
# output result and statistics
open OUT, ">$ARGV[1].cutpoints" or die( "$!" );
print OUT $result;
close OUT;

open OUT, ">$ARGV[1].cutpoints.stat" or die( "$!" );
print OUT "Alignable\t$alignable\nDiscarded\t$discarded\nPass\t$pass\n";
close OUT;

sub analyze_segment {
	my $rid = shift;
	my $records = shift;
	return unless $rid;		## empty record
	++ $alignable;

	## extract R1 and R2
	my (@R1s, @R2s);
	foreach my $hit ( @$records ) {
		if( $hit->[0] & 64 ) {
			push @R1s, $hit;
		} elsif( $hit->[0] & 128 ) {
			push @R2s, $hit;
		} else {	 ## unknown flag, discard it
			print STDERR "E0\tunacceptable flag ($hit->[0]) for $rid!\n";
		}
	}
	if( $#R1s < 0 && $#R2s < 0 ) {	## no valid hits
		++ $discarded;
		return;
	}

	## check segment numbers, I can rescue 1 end, which is usually caused due to too-short of the other end
	my $end_to_rescue = '';
	if( $#R1s < 0 || $#R2s < 0 ) {
		if( $#R1s==0 && $rescueMode ) {
			$end_to_rescue = $R1s[0];
		} elsif( $#R2s==0 && $rescueMode ) {
			$end_to_rescue = $R2s[0];
		} else {
			++ $discarded;
			return;
		}
		my ($chr, $cutpoint, $strand) = extract_cutpoint( $end_to_rescue );
		if( $chr eq '' ) {
			++ $discarded;
			return;
		}
		## check duplicate
#		my $extra = ($rid=~/#(\S+)/) ? $1 : '';
#		my $key = "$chr:$cutpoint:$strand:$extra";
#		if( exists $met{$key} ) {
#			++ $duplicate;
#		} else {
#			$met{$key} = 1;
			++ $pass;
			$result .= "Cutpoint\t$chr\t$cutpoint\t$strand\t$rid\n";
#		}
		return;
	}

	## both R1 and R2 have hits
	if( $#R1s==0 && $#R2s==0 ) {
		my ($chr1, $cutpoint1, $strand1) = extract_cutpoint( $R1s[0] );
		my ($chr2, $cutpoint2, $strand2) = extract_cutpoint( $R2s[0] );
		if( $chr1 eq '' && $chr2 eq '' ) {
			++ $discarded;
			return;
		}

		my ($key1, $key2, $code);
		if( $chr1 eq '' ) {
			$key1 = "$chr2\t$cutpoint2\t$strand2";
			$code = "Cutpoint";
		} elsif( $chr2 eq '' ) {
			$key1 = "$chr1\t$cutpoint1\t$strand1";
			$code = "Cutpoint";
		} else {
			if( $chr1 lt $chr2 ) {
				$key1 = "$chr1\t$cutpoint1\t$strand1\t$chr2\t$cutpoint2\t$strand2";
			} elsif( $chr1 eq $chr2 && $cutpoint1<$cutpoint2 ) {
				$key1 = "$chr1\t$cutpoint1\t$strand1\t$chr2\t$cutpoint2\t$strand2";
			} else {
				$key1 = "$chr2\t$cutpoint2\t$strand2\t$chr1\t$cutpoint1\t$strand1";
			}
			$code = "Pair";
		}

		if( $code eq 'Cutpoint' ) {
#			if( exists $met{$key1} ) {
#				++ $duplicate;
#				return;
#			} else {
#				$met{$key1} = 1;
				++ $pass;
#			}
		} elsif( $code eq 'Pair' ) {
#			if( exists $met{$key1} || exists $met{$key2} ) {
#				++ $duplicate;
#				return;
#			} else {
#				$met{$key1} = 1;
				++ $pass;
#			}
		}
		$result .= "$code\t$key1\t$rid\n";
	} else {	## something is wrong, too many segments
		++ $discarded;
		return;
	}
}

