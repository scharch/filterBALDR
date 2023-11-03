#!/usr/bin/env perl

=head 1 SUMMARY

 filterBALDR.pl
 Usage: filterBALDR.pl sample_sheet.csv [cellregex]

 This script filters the output of BALDR single cell V(D)J assembly to remove
     short contigs and those with low read counts (5-fold drop-off) and outputs
     the remaining contigs in fasta format for use with SONAR, IMGT, etc.
 Assumes the presence of BALDR output for both heavy and light chains, with the
     default names. Output is `filtered.fa` in the working directory.

 Parameters:
     sample_sheet.csv     Sample sheet used as input to BALDR. Used to convert
                              index (column 4) to cell id (column 7). Assumes
                              that the first column of the BALDR summary includes
                              the path "<flow_cell>/<lane>".
     cellregex            Regular expression pattern for extracting the index
                              from the BALDR file name.
                              [default: "UDP\d{4}_UDP\d{4}"]

 Created by Chaim A. Schramm 2019-08-28
 Edited input options to fit better with standard usage by CA Schramm, 2021-06-23.
 Changed to filtering read counts by plate, instead of by well, CAS 2021-07-12.
 Changed how <sample_sheet> is processed to allow handling of multiple
     flowcells/lanes at once from the top-level directory by CA Schramm 2021-08-11.
 Changed threshold to 1% of median (instead of 10%) by CA Schramm 2021-08-11.
 Fixed plate identification by CA Schramm 2023-11-02.

 Copyright (c) 2019-2023 Vaccine Research Center, NIAID, National Institutes of Health, USA.
 All rights reserved.
=cut

use strict;
use diagnostics;
use Pod::Usage;
use Statistics::Basic qw/median/;

if ($#ARGV<0 || $ARGV[0] =~ /-h/) { pod2usage(1); }

my ($sampleSheet, $cellPattern) = @ARGV;
$cellPattern = 'UDP\d{4}_UDP\d{4}' unless -defined($cellPattern);

#load sample sheet
my %lookUp;
open SAMPLE, $sampleSheet or die "Can't read from $sampleSheet: $!\n";
while(<SAMPLE>) {
	chomp;
	my @arr = split(/,/);
	$lookUp{ $arr[2] }{ $arr[3] }{ $arr[4] } = $arr[7];
}
close SAMPLE;


my %contigs;
for my $file ("Results_IGH_rank_all.txt", "Results_IGKL_rank_all.txt") {

	my %wells = ();
	my %rank1 = ();

	open IN, $file or die "Can't read from $file: $!\n";

	my $header = <IN>; #discard
	while (<IN>) {

		next if /^$/; #skip random blank lines
		chomp;
		my @row = split(/\t/);
		next if $#row < 5; #skip partial lines from cells with no data
		next if $row[12] eq "IGHV7-40*03"; #ignore bogus contigs from IgBLAST 1.15

		#some prefiltering for junk contigs
		next unless $row[48] =~ /^(IG[KL]J|J[KL])/ || $row[49] =~ /^(IGHJ|JH)/; #Top J gene match; column wobbles depending on presence or absence of D
		next unless $row[13] >= 50; #V alignment length
		next unless length( $row[-1] ) >= 75; #V(D)J sequence
		next if $row[10] eq "-"; #no supporting reads???

		my ($flowcell, $lane, $index) = $row[0] =~ /(\d{6}_(?:M|A|VH)\d+_\d+_(?:000000000-)?[^\/]+)\/(\d)\/.+($cellPattern)/;
		my $cellID = "$flowcell-$lane-$index";
  		my $plate  = "overflow"; #will threshold all data with no plate identification together

		if (exists $lookUp{$flowcell}{$lane}{$index}) {
			$cellID = $lookUp{$flowcell}{$lane}{ $index };
   			($plate) = $cellID =~ /^(\d+)\-?[A-Pa-p]\d+$/;
      			if($plate eq "") {
	 			warn("Can't identify plate for cell $cellID, will threshold separately...");
     				$plate = "overflow";
	 		}
		} else {
			warn("Couldn't find $flowcell-$lane-$index in $sampleSheet, using as is...");
		}

		$row[0] = $cellID; #to avoid recalculating

		push @{$wells{$plate}}, \@row;
		push @{$rank1{$plate}}, $row[10] if $row[3] eq "1";

	}

	close IN;


	#go back through the plates and filter by read support
	for my $plate (sort keys %rank1){
		my $med = median(@{$rank1{$plate}});
		for my $line (@{$wells{$plate}}) {
			push @{$contigs{$line->[0]}}, $line if $line->[10] >= 0.01 * $med;
		}
	}
}

open OUT, ">filtered.fa" or die "Can't write to filtered.fa: $!\n";
for my $cell (sort keys %contigs) {
	for ( my $i=0; $i<scalar(@{$contigs{$cell}}); $i++ ) {
		my $num = $i+1;
		print OUT ">$cell.$num cell_id=$cell duplicate_count=$contigs{$cell}[$i][10]\n$contigs{$cell}[$i][-3]\n"; #10 is Bowtie2_VDJ count; -3 is raw full contig field
	}
}
close OUT;
