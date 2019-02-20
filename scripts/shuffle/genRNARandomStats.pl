#!/usr/bin/perl -I/tmp/ViennaRNA-1.5/Perl/blib/arch -I/tmp/ViennaRNA-1.5/Perl/blib/lib 
#-I/home/stanley/ViennaRNA-1.5/Perl/blib/arch -I/home/stanley/ViennaRNA-1.5/Perl/blib/lib -I/home/stanley/Shuffle-1.4/blib/lib -I/home/stanley/Statistics-Basic-0.41.3/blib/lib
#-I/home/stanley/ViennaRNA-1.5/Perl/blib/arch -I/home/stanley/ViennaRNA-1.5/Perl/blib/lib -I/home/stanley/Shuffle-1.4/blib/lib -I/home/stanley/Statistics-Descriptive-2.6/blib/lib

############################################################################
# AUTHOR:  	Stanley NG Kwang Loong, stanley@bii.a-star.edu.sg
# DATE:		31/07/2005
# VERSION:	1.0
############################################################################
use Algorithm::Numerical::Shuffle qw /shuffle/;
use Statistics::Basic::StdDev;
use Statistics::Basic::Mean;
use RNA;
$ENV{UNBIAS}=1;

use POSIX qw(log);
use strict;
srand(time|$$);


############################################################################
# Global Parameters and initialization if any.
############################################################################

my $mfeFile="";
my $inFile="&STDIN";
my $outFile="&STDOUT";
my $numRandomSeqs; #number of random sequences per set 
my $numSeqs; 
my $usage = "USAGE:	perl genRNARandomStats.pl -n \"number\" -i <input file> -o <output file> -m <mfe file>\n";
my $statsID=-1;
my @templines;
my @amfe;
my @aseq;
my @astruct;
my @aQ;
my @aD;
my @aBP;
my @aSS;
my @cmfe;
my @cseq;
my @cstruct;
my @cQ;
my @cD;
my @cBP;
my @cSS;

############################################################################
# File IO
# Parse the command line.
############################################################################

foreach my $a (0..$#ARGV) {

	if ($ARGV[$a] eq "-n") {	
		$numRandomSeqs = $ARGV[$a+1];
	}        
	elsif ($ARGV[$a] eq "-i") {
		$inFile = ($ARGV[$a+1] eq "") ? "&STDIN" : $ARGV[$a+1];	
	}
	elsif ($ARGV[$a] eq "-o") {	
		$outFile = ($ARGV[$a+1] eq "") ? "&STDOUT" : $ARGV[$a+1];
		printf(STDERR "$outFile\n");
		if ( -e $outFile) {
			open (OUTFILE, "<$outFile") or die( "Cannot open input file $outFile: $!" );
			while (my $line = <OUTFILE>) {
				chomp($line);
				if ($line =~ m/^[0-9]/) {
					push(@templines, $line);
					$statsID++;
				}
			}
			close (OUTFILE);
		}
	}
	elsif ($ARGV[$a] eq "-m") {	
		$mfeFile = $ARGV[$a+1];
	}
	else { }
}

if(!defined($inFile) || !defined($outFile) || !defined($mfeFile) || !defined($numRandomSeqs)) { 
	die ("$usage\n"); 
}
open (INFILE, "<$inFile") or die( "Cannot open input file $inFile: $!" );
open (MFEFILE, "<$mfeFile") or die( "Cannot open input file $mfeFile: $!" );

# ID Mean SD
$numSeqs = 0;
while (my $line = <MFEFILE>) {
	chomp($line);
	if ($line =~ m/^>/) { }
	
	elsif ($line =~ m/^[AaCcUuGg]/) {
		$aseq[$numSeqs] = $line;
	}
	# Fasta Third Line i.e. RNA secondary structure and MFE
    elsif ($line =~ m/\)$/) {
		$line =~ s/\( /\(/;
		($astruct[$numSeqs], $amfe[$numSeqs]) = split(/ /, $line);
		$amfe[$numSeqs] =~ s/[()]//g;
		($aQ[$numSeqs], $aD[$numSeqs], $aBP[$numSeqs], $aSS[$numSeqs]) = rnaAnalysis($aseq[$numSeqs], length($aseq[$numSeqs]), $astruct[$numSeqs], $amfe[$numSeqs]);
		$numSeqs++;
		
	}			
    else {}
  
}#end of while loop
close (MFEFILE) or die( "Cannot close input file $mfeFile: $!" );


open (OUTFILE, ">$outFile") or die ("Cannot open output file $outFile: $!");
print (OUTFILE "ID");
print (OUTFILE "\tX\tRMean\tRSD\tZ\tP-MFE\tP-Z");
print (OUTFILE "\tX\tRMean\tRSD\tZ\tP-Q\tP-Z");
print (OUTFILE "\tX\tRMean\tRSD\tZ\tP-D\tP-Z");
print (OUTFILE "\tX\tRMean\tRSD\tZ\tP-PB\tP-Z");
print (OUTFILE "\tX\tRMean\tRSD\tZ\tP-TD\tP-Z");
map { print (OUTFILE "\n$_"); } @templines; 

my $i = 0; #index for $numRandomSeqs
$numSeqs = 0;
while (my $line = <INFILE>) { # Read line by line.
	chomp($line);
	if ($line =~ m/^>/) { }
	
	elsif ($line =~ m/^[AaCcUuGg]/) {
		 
		$cseq[$i] = $line if ($numSeqs > $statsID);
	}
	
	# Fasta Third Line i.e. RNA secondary structure and MFE
    elsif ($line =~ m/^[.(]/) {
    	if ($numSeqs > $statsID) {
			$line =~ s/\( /\(/;
			($cstruct[$i], $cmfe[$i]) = split(/ /, $line);
			$cmfe[$i] =~ s/[()]//g;

			($cQ[$i], $cD[$i], $cBP[$i], $cSS[$i]) = rnaAnalysis($cseq[$i], length($cseq[$i]), $cstruct[$i], $cmfe[$i]);
		}
		$i++;
		
		#Finished collecting $numRandomSeqs sequences
		if ($i == $numRandomSeqs) {
			if ($numSeqs > $statsID) {
				printf (OUTFILE "\n%u", $numSeqs+1);
				my ($mean, $sd, $z, $p, $pz);
				($mean, $sd, $z, $p, $pz) = computeStats(\@cmfe, $amfe[$numSeqs]);
				printf (OUTFILE "\t%.2f\t%.4f\t%.4f\t%s\t%.4f\t%s", $amfe[$numSeqs], $mean, $sd, $z, $p, $pz);
				($mean, $sd, $z, $p, $pz) = computeStats(\@cQ, $aQ[$numSeqs]);
				printf (OUTFILE "\t%.2f\t%.4f\t%.4f\t%s\t%.4f\t%s", $aQ[$numSeqs], $mean, $sd, $z, $p, $pz);
				($mean, $sd, $z, $p, $pz) = computeStats(\@cD, $aD[$numSeqs]);
				printf (OUTFILE "\t%.2f\t%.4f\t%.4f\t%s\t%.4f\t%s", $aD[$numSeqs], $mean, $sd, $z, $p, $pz);
				($mean, $sd, $z, $p, $pz) = computeStats(\@cBP, $aBP[$numSeqs]);
				printf (OUTFILE "\t%.2f\t%.4f\t%.4f\t%s\t%.4f\t%s", $aBP[$numSeqs], $mean, $sd, $z, $p, $pz);

				my @cstructtree = map { RNA::make_tree(RNA::expand_Full($_)) } @cstruct;    
				my $astructtree = RNA::make_tree(RNA::expand_Full($astruct[$numSeqs]));
				my @aedit_distance = map { RNA::tree_edit_distance($astructtree, $_) } @cstructtree;
				my $aedit_distance_mean = new Statistics::Basic::Mean(\@aedit_distance)->query;

				@cstructtree = shuffle(@cstructtree);
				my $chosentree = pop(@cstructtree);
				my @cedit_distance = map { RNA::tree_edit_distance($chosentree, $_) } @cstructtree;

				($mean, $sd, $z, $p) = computeTreeStats(\@cedit_distance, $aedit_distance_mean);
				$pz = 0.0;
				printf (OUTFILE "\t%.2f\t%.4f\t%.4f\t%s\t%.4f\t%s", $aedit_distance_mean, $mean, $sd, $z, $p, $pz);
			}
			$i = 0;	
			$numSeqs++;			
			last if (scalar(@amfe) == $numSeqs);
		}
	}			
    else { }
  
}#end of while loop

print(OUTFILE "\n");

close (INFILE) or die( "Cannot close input file $inFile: $!" );
close (OUTFILE) or die( "Cannot close output file $outFile: $!");

sub pzCount {
	my ($nativeZ, $times, $data) = @_;
	my $size = 1000;
	
	return 0 if (scalar(@$data) < $size+1);
	
	my $R = 0;
	my $times_size = sprintf("%d", (scalar(@$data) - 1)/$size);
	
	while (1) {
	
		my @shuffle = shuffle(@$data);
		my $first = 0;
		my $last = $size - 1;
		my $index = scalar(@$data) - 1;
		foreach my $i(0..$times_size-1) {
			my @short = @shuffle[$first..$last];
			my $mean = new Statistics::Basic::Mean(\@short)->query;
			my $sd = new Statistics::Basic::StdDev(\@short)->query; 
			$R++ if (($shuffle[$index] - $mean) < ($nativeZ * $sd));
			return $R if (--$times == 1);
			$first += $size;
			$last += $size;
			$index -= 1;
		}
	}

}

sub computeStats {

	my ($data, $value) = @_;
	my $mean = new Statistics::Basic::Mean($data)->query;
	my $sd = new Statistics::Basic::StdDev($data)->query;
	my $p = 0;
	map { $p++ if ($_ < $value); } @$data;
	my $z = "undef";
	my $pz = "undef";
	if ($sd > 0) {
		$z = sprintf("%.4f", ($value - $mean)/$sd);
		$pz = 0.0;
	}
	
	return ($mean, $sd, $z, $p/scalar(@$data), $pz);
}

sub pzTreeCount {
	my ($nativeZ, $times, $data) = @_;
	my $size = 1000;
	
	return 0 if (scalar(@$data) < $size+1);
	
	my $R = 0;
	my $times_size = sprintf("%d", (scalar(@$data) - 1)/$size);
	
	while (1) {
	
		my @shuffle = shuffle(@$data);
		my $first = 0;
		my $last = $size - 1;
		my $index = scalar(@$data) - 1;
		foreach my $i(0..$times_size-1) {
			my @sample = @shuffle[$first..$last];

			my @aedit_distance = map { RNA::tree_edit_distance($shuffle[$index], $_) } @sample;
			my $aedit_distance_mean = new Statistics::Basic::Mean(\@aedit_distance)->query;
			
			my $chosentree = pop(@sample);
			my @cedit_distance = map { RNA::tree_edit_distance($chosentree, $_) } @sample;			
			my $mean = new Statistics::Basic::Mean(\@cedit_distance)->query;
			my $sd = new Statistics::Basic::StdDev(\@cedit_distance)->query; 

			$R++ if (($aedit_distance_mean - $mean) < ($nativeZ * $sd));
			return $R if (--$times == 1);
			$first += $size;
			$last += $size;
			$index -= 1;
		}
	}
}

sub computeTreeStats {

	my ($data, $value) = @_;

	my $mean = new Statistics::Basic::Mean($data)->query;
	my $sd = new Statistics::Basic::StdDev($data)->query;
	my $p = 0;
	map { $p++ if ($_ < $value); } @$data;
	my $z = "undef";
	if ($sd > 0) {
		$z = sprintf("%.4f", ($value - $mean)/$sd);
	}
	
	return ($mean, $sd, $z, $p/scalar(@$data));
}


sub pairwiseTreeDistance {
	my ($data, $value) = @_;
    my @datatree = map { RNA::make_tree(RNA::expand_Full($_)) } @$data;    
    my $valuetree = RNA::make_tree(RNA::expand_Full($value));
    
	my @temp = map { RNA::tree_edit_distance($valuetree, $_) } @datatree;
	my $tempmean = new Statistics::Basic::Mean(\@temp)->query;

	@datatree = shuffle(@datatree);
	my $chosentree = pop(@datatree);
	@temp = ();
	@temp = map { RNA::tree_edit_distance($chosentree, $_) } @datatree;
	
	return(computeStats(\@temp, $tempmean));
}

sub rnaAnalysis {
	my ($seq, $seqLen, $struct, $mfe) = @_; 
	my $bp = $struct =~ tr/(//; 
	my $Q = 0;
	my $D = 0;
	my $VI = 0;

	$RNA::pf_scale = exp((-1)*1.2*$mfe/(0.6163207755*$seqLen)) if ($seqLen > 1000);

	# compute partition function and pair pobabilities matrix
	RNA::pf_fold($seq);   				
	# compute sum-of-entropy and bp-distance
	foreach my $j (1..$seqLen-1) {
		foreach my $k ($j+1..$seqLen) {
			my $p = RNA::get_pr($j, $k); # points to the computed pair probabilities
			if ($p > 0) {
				$Q += $p*log($p);
				$D += $p*(1 - $p);
			}
		}
	}
	$Q *= -1/log(2);
	return ($Q, $D, $bp, 0);

}

exit;

