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
	print STDERR "\nUsage: $0 <border.sam> <out.prefix>";
	print STDERR "\nExtract duplcate-removed ends. The Input sam MUST be sorted by name.\n\n";
	exit 2;
}

my ($alignable, $discarded, $pass) = (0,0,0);
my $result = '';

## start analysis
my $lastID = '';
my @lastRead = ();
my $cnt = 0;
if( $ARGV[0] =~ /bam$/ ) {
	open IN, "samtools view $ARGV[0] |" or die( "$!" );
} else {
	open IN, "$ARGV[0]" or die( "$!" );
}
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

	if( $#$records == 0 ) {	## 1 segments
		++ $alignable;

		my ($chr, $cutpoint, $strand) = extract_cutpoint( $records->[0] );
		if( $chr eq '' ) {
			++ $discarded;
			return;
		}
		## check duplicates
#		my $extra = ($rid=~/#(\S+)/) ? $1 : '';
#		my $key = "$chr:$cutpoint:$strand:$extra";
#		print STDERR "$rid\t$key\n";
#		if( exists $met{$key} ) {
#			++ $duplicate;
#		} else {
#			$met{$key} = 1;
			++ $pass;
			$result .= "Cutpoint\t$chr\t$cutpoint\t$strand\t$rid\n";
#		}
	} else {	## more than 1 segments, discard
		++ $alignable;
		++ $discarded;
	}
}

