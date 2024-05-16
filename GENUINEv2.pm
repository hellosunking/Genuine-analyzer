#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

package GENUINEv2;
use FindBin qw/$Bin/;

use Exporter 'import';
@EXPORT_OK = qw/get_clip_size_from_cigar extract_cutpoint $spikeSeq $spikeSize $min_mapQ $min_mapped_ratio $max_unmapped_allowed SmithWaterman visualizeSW revcmp loadGenome/;

## Spike-in information
our $spikeSeq  = 'ACGGTGGACCGATGATCATCGGTCCACCGT';
our $spikeSize = length( $spikeSeq );

## parameters that is NOT configurable by the users
our $min_mapQ = 30;
our $min_mapped_ratio = 0.75;	## the read must be covered at least 90% to consider
our $min_clip_size    = 10;		## minimum bp to be considered as a clip; STAR sometimes give XXXM1S thing
our $max_unmapped_allowed = 20;	## minimum unmapped bp to discard the alignment, this is to rescue 100M15S thing
our $maxDist = 1000;			## for check_pair

sub get_clip_size_from_cigar {
	my $cigar = shift;

	my ($leftClip, $rightClip, $size) = ( 0, 0, 0 );
	if( $cigar =~ s/^(\d+)[HS]// ) {
		$leftClip = $1;
	}
	if( $cigar =~ s/(\d+)[HS]$// ) {
		$rightClip = $1;
	}

	$cigar =~ s/([HSMIDN])/$1:/g;
	my @seg = split /:/, $cigar;
	foreach my $c ( @seg ) {
		if( $c =~ /^(\d+)([HSMIDN])$/ ) {
			if ( $2 eq 'M' || $2 eq 'D' ) {	## match/deletion, count it
				$size += $1;
			} elsif( $2 eq 'I' ) {	## insertion, skip
				next;
			} else {
				print STDERR "ERROR 255: invalid cigar ($cigar) for $rid!\n";
				return (0,0,0);
			}
		} else {
			print STDERR "ERROR 255: invalid cigar ($cigar) for $rid!\n";
			return (0,0,0);
		}
	}

	return ($leftClip, $rightClip, $size);
}

sub extract_cutpoint {
	my $hit = shift;

	my $chr   = $hit->[1];
	my $start = $hit->[2];
	my $cigar = $hit->[-1];
	my $strand = ($hit->[0] & 16 ) ? '-' : '+';
	my $mapQ = $hit->[3];
	my ($leftClip, $rightClip, $templateSize) = get_clip_size_from_cigar( $cigar );
	if( $templateSize ==0 ) {       ## error occurs
		++ $discarded;
		return ('', '', '');
	}## check read integrity
	if( $templateSize < ($templateSize+$leftClip+$rightClip) * $min_mapped_ratio ) {
		++ $discarded;
		return ('', '', '');
	}
	my $cutpoint;## check duplicates
	if( $strand eq '+' ) {
		$cutpoint = $start;
	} else {
		$cutpoint = $start + $templateSize - 1;
	}
	return ($chr, $cutpoint, $strand);
}

#local alignment by using Smith-Waterman algorithm
sub SmithWaterman {
    #get the options
    my ($seq1, $seq2) = @_;

	#set the default value of key paramters
	my $match   =  1;
	my $unmatch = -1;
	my $space   = -2;

    #get the dimension of matrix
    my $row_num = (length $seq1) + 2;
    my $col_num = (length $seq2) + 2;

    #create 2 matrices, one for score and one for arrow
    my (@matrix_score, @matrix_arrow);
    push @matrix_score, [ (0) x $col_num ] foreach 1 .. $row_num;
    push @matrix_arrow, [ (0) x $col_num ] foreach 1 .. $row_num;

    #step1----------------------------------------------------------------------
    #for score and arrow matrix, fill the seq of row 1 and col 1, and fill the score/arrow of row 2 and col 2
    #note : in arrow matrix, -1 mean left, 1 mean up, 0 mean upper left, and 2 mean empty
    foreach my $i (0 , 1) {
        foreach my $j (2 .. $#{$matrix_score[$i]}) {
            $matrix_score[$i][$j] = $i == 0 ? substr $seq2, $j - 2, 1 : 0;
            $matrix_arrow[$i][$j] = $i == 0 ? substr $seq2, $j - 2, 1 : 2;
        }
    }
    foreach my $j (0 , 1) {
        foreach my $i (2 .. $#matrix_score) {
            $matrix_score[$i][$j] = $j == 0 ? substr $seq1, $i - 2, 1 : 0;
            $matrix_arrow[$i][$j] = $j == 0 ? substr $seq1, $i - 2, 1 : 2;
        }
    }

    #step2----------------------------------------------------------------------
    #for score and arrow matrix, fill the other cells, respectively
    foreach my $i (2 .. $#matrix_score) {
        foreach my $j (2 .. $#{$matrix_score[$i]}) {
            my $road_1_score = $matrix_score[$i][$j-1] + $space; #reach each cell from the left
            my $road_2_score = $matrix_score[$i-1][$j] + $space; #reach each cell from the up
            my $road_3_score = $matrix_score[$i][0] eq $matrix_score[0][$j] ?
                               $matrix_score[$i-1][$j-1] + $match : $matrix_score[$i-1][$j-1] + $unmatch; #reach each cell from the upper left
            if ($road_1_score > $road_2_score) {
                $matrix_score[$i][$j] = $road_1_score > $road_3_score ? $road_1_score : $road_3_score;
            }
            else{
                $matrix_score[$i][$j] = $road_2_score > $road_3_score ? $road_2_score : $road_3_score;
            }
            #if there are multiple paths for the maximum score, road 3 will be used first, road 2 will be used second, and road 1 will be used last
            $matrix_arrow[$i][$j] = $matrix_score[$i][$j] == $road_3_score ? 0 : ($matrix_score[$i][$j] == $road_2_score ? 1 : -1);
            #if the cell of the score matrix is negative, replace it with 0, and the corresponding arrow matrix cell is empty
            ($matrix_score[$i][$j], $matrix_arrow[$i][$j]) = (0, 2) if $matrix_score[$i][$j] < 0;
        }
    }

    #step3----------------------------------------------------------------------
    #get a maximum alignment score
    my ($row, $col) = (1, 1);
    foreach my $i (2 .. $#matrix_score) {
        foreach my $j (2 .. $#{$matrix_score[$i]}) {
            ($row, $col) = ($i, $j) if $matrix_score[$i][$j] > $matrix_score[$row][$col];
        }
    }
    my $max_align_score = $matrix_score[$row][$col];

    #get two aligned sequences
    my ($align_seq_1, $align_seq_2) = ('', '');
    my ($i, $j) = ($row, $col);
    while ($matrix_score[$i][$j]) {
        if ($matrix_arrow[$i][$j] == -1) {
            $align_seq_1 .= '-';
            $align_seq_2 .= $matrix_arrow[0][$j];
            $j--;
        }
        elsif ($matrix_arrow[$i][$j] == 1) {
            $align_seq_1 .= $matrix_arrow[$i][0];
            $align_seq_2 .= '-';
            $i--;
        }
        else{
            $align_seq_1 .= $matrix_arrow[$i][0];
            $align_seq_2 .= $matrix_arrow[0][$j];
            $i--;
            $j--;
        }
    }
    $align_seq_1 = reverse $align_seq_1;
    $align_seq_2 = reverse $align_seq_2;

	# calculate the aligned sequence
	-- $row; -- $col;
	-- $i; -- $j;
	my $ref = substr($seq1, 0, $i);
	$ref .= $align_seq_1;
	$ref .= substr($seq1, $row);
	my $qry = ' ' x ($i-$j);
	if( $j ){ $qry .= substr($seq2, 0, $j); }
	$qry .= $align_seq_2;
	if( $col != length($seq2) ) {
		$qry .= substr($seq2, $col);
	}
#	print STDERR "ref: $ref\nqry: $qry\n";

	## NOTE: I will recalculate the score here because the indels are scored as -2 during searching
	## and I will also ignore the mismatch to NGG
	my $mismatch = 0;
	foreach my $sub ( 0 .. length($qry)-1 ) {
		my $qq = substr($qry, $sub, 1);
		next if $qq eq ' ' || $qq eq 'N';	## for the NGG suffix in Cas9
		$mismatch += $qq ne substr($ref, $sub, 1);
	}
#	print STDERR "newScore=$newScore, OLD Score=$max_align_score\n";

    #return the result
    return ($mismatch, $max_align_score, $align_seq_1, $align_seq_2, $ref, $qry);
}

#step2 : add align line between the two sequences
sub visualizeSW {
	my $ref = shift;
	my $qry = shift;	## design
	my ($match, $mismatch, $insertion, $deletion) = (0,0,0,0);

	#step2 : add align line between the two sequences
	my $maxima = (length($qry) < length($ref)) ? length($qry) : length($ref);
	my $s;
	for($s=0; $s<$maxima; ++$s) {
		last if substr($qry, $s, 1) ne " ";
	}
	if( $s > $maxima ) {
		print STDERR "Something is wrong: ref=$ref, qry=$qry\n";
		return 'NO_ALIGN';
	}

	my @align_line;
	for(my $i=$s; $i<$maxima; ) {
		my $a = substr($ref, $i, 1);
		my $q = substr($qry, $i, 1);

		if( $q ne "-" ) {	## non-insertion
			if($a eq $q) {
				push @align_line, '.';
				++ $match;
			} else {
				push @align_line, $a;
				if( $a eq '-' ) {
					++ $deletion;
				} else {
					++ $mismatch;
				}
			}
			++ $i;
		} else {	## insertion, tough
			my $insert = $a;
			while( $i<$maxima ) {
				++ $i;
				$q = substr($qry, $i, 1);
				if( $q ne "-" ) {	## end of insertion
					last;
				} else {
					$insert .= substr($ref, $i, 1);
				}
			}
			$align_line[-1] .= $insert;
##			++ $insertion;	## consider nII as 1 insertion
			$insertion += length($insert);	## consider nII as 2 insertions
		}
	}

	return (join(";", @align_line), $match, $mismatch, $insertion, $deletion);
}

sub revcmp {
	my $s = shift;
	$s = reverse $s;
	$s =~ tr/ACGT/TGCA/;
#	$s =~ tr/RYMKHBVD/YRKMDVBH/;
	return $s;
}

sub loadGenome {
	my $gid = shift;
	my $file = ($gid=~/\.fa$/) ? $gid : "$Bin/genome/$gid.fa";

	my %g;
	print STDERR "Loading genome: $gid ($file)\n";
	open GENOME, "$file" or die("$!");
	my $chr = 'NA';
	while( <GENOME> ) {
		chomp;
		if( /^>(\S+)/ ) {
			$chr = $1;
		} else {
			$g{$chr} .= uc $_;
		}
	}
	close GENOME;
	return \%g;
}

1;

