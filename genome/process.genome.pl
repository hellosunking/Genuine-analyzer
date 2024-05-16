#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.fa[.gz]> [in.fa ...]\n\n";
	exit 2;
}

foreach my $file ( @ARGV ) {
	if( $file =~ /\.gz$/ ) {
		open IN, "zcat $file |" or die( "$!" );
	} else {
		open IN, "$file" or die( "$!" );
	}

	my $badchr = 1;
	while( <IN> ) {
		if( /^>(\S+)/ ) {
			$badchr = ($1=~/_/);
		}
		print unless $badchr;
	}
	close IN;
}

