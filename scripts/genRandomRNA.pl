#!/usr/bin/perl 
#-w -I/home/stanley/bioperl-1.4/blib/lib -I/home/stanley/Shuffle-1.4/blib/lib/  

############################################################################s
# AUTHOR:  	Stanley NG Kwang Loong, stanley@bii.a-star.edu.sg
# DATE:		31/07/2005
# FILENAME: gen_random_rna.pl
# VERSION:	1.0
# DESCRIPTION: Generate some (pseudo-)random rna sequences
#				
# USAGE:	perl genRandomRNA.pl -n 'num repeats' -m 'method' < <input file> > <output file> 
#			n: Number of repeats per sequence. For example 10 will produce 10 shuffled sequences per original sequence.
#			m: Method of shuffling (m: mononucleotid, d: dinucleotide z: zero-order markov, f:first-order markov ).
############################################################################
use strict;
use Getopt::Std;
use Bio::SeqIO;
use Algorithm::Numerical::Shuffle qw /shuffle/;

my %gl_monomers = ('A' => 0,'C' => 0, 'G' => 0, 'U' => 0);
my %gl_dimers = ('AA' => 0, 'AC' => 0, 'AG' => 0, 'AU' => 0,
	                 'CA' => 0, 'CC' => 0, 'CG' => 0, 'CU' => 0,
	                 'GA' => 0, 'GC' => 0, 'GG' => 0, 'GU' => 0,
                 	 'UA' => 0, 'UC' => 0, 'UG' => 0, 'UU' => 0);

# Seed the random number generator.
# time|$$ combines the current time with the current process id
srand(time|$$);

my $usage = "USAGE:	perl genRandomRNA.pl -n 'num repeats' -m 'method' < <input file> > <output file>\n";

if (scalar(@ARGV) < 2) {
  print(STDERR $usage);
  exit(1);
}

# set options
my %opt=();
getopts("n:m:",\%opt);

my $repeat = $opt{"n"};
my $method = $opt{"m"};

if ($method !~ /[m|d|z|f]/ || $repeat !~ /\d+/) { die ("$usage\n"); }

# Get sequences from fasta file
my $in = new Bio::SeqIO(-fh => \*STDIN, -format => "fasta");

# Process all sequences, upper case and convert to ACGU
while (my $obj = $in->next_seq) {

	my $id = $obj->id;
	my $seq = uc($obj->seq);
	$seq =~ s/T/U/g;
	
	# Mono-nucleotide shuffle (permute) sequence generation
	if ($method eq "m") {
		for my $i (0..$repeat-1) {
			my @tmp = shuffle(split (//,$seq));
			$seq = join('', @tmp);
			print (">$id-$i\n$seq\n");
		}
	}
	elsif ($method eq "d") {
		for my $i (0..$repeat-1) {
			my $seq = altschulEriksonDinuclShuffle($seq);
			print (">$id-$i\n$seq\n");
		}
	}
	elsif ($method eq "z") {
		for my $i (0..$repeat-1) {
			my $seq = genZMMRNASet($seq);
			print (">$id-$i\n$seq\n");
		}
	}
	elsif ($method eq "f") {
		for my $i (0..$repeat-1) {
			my $seq = genFMMRNASet($seq);
			print (">$id-$i\n$seq\n");
		}
	}
	else { }
}

############################################################################
# Di-nucleotide shuffle (permute) sequence generation
############################################################################
sub altschulEriksonDinuclShuffle{

	my ($seq) = @_;
	my $lastIndex = length($seq)-1;
	my $lastCh = substr($seq, $lastIndex, 1);
	my %edgeList;

	#Create Edge Lists
	for my $i (0..$lastIndex-1) {
		$edgeList{substr($seq, $i, 1)} .= substr($seq, $i+1, 1);
	}

	# Choose three edges that form a subgraph ending with $lastCh
	NOTFOUND: while (1) {
		my @edgeStr;
		my %D;	
		
		#Shuffle Edge Lists that doesn't start with $lastCh, and pick their last item for creating edges
		foreach my $a (keys (%edgeList)) {	
			if ($a ne $lastCh) {
				my @temp = shuffle(split(//, $edgeList{$a})); 
				$edgeList{$a} = join('', @temp);
				push(@edgeStr, $a);
				my $b = substr($edgeList{$a}, @temp-1, 1);
				push(@edgeStr, $b);
				$D{$a} = ($b eq $lastCh) ? 1 : 0;
			}		
		}
		
		#Does the chosen edges individually end with $lastCh ?
		my $sum;
		foreach my $DValue (values (%D)) {
			$sum += $DValue;
		}
		next NOTFOUND if ($sum == 0);
		
		#Try to form a subgraph ending with $lastCh
		for my $i (0..keys(%D)) { 
			for (my $j = 0; $j < @edgeStr; $j+=2) {
				$D{$edgeStr[$j]} = 1 if (($edgeStr[$j+1] eq $lastCh) || ($D{$edgeStr[$j+1]} == 1));
			}
		}
		
		#Does the chosen edges form a subgraph ending with $lastCh ?
		foreach my $i (keys (%D)) {	
			 next NOTFOUND if (($i ne $lastCh) && ($D{$i} == 0));		
		}

		last;
	}
	
	#Construct the eulerian path
	my $shuf_seq = substr($seq, 0, 1);	
	my $i = length($seq);
	
	$edgeList{$lastCh} = join('', shuffle(split(//, $edgeList{$lastCh}))) if (defined $edgeList{$lastCh});
				
	while (--$i) {
		my $prevCh = substr($shuf_seq, length($shuf_seq) - 1, 1);
		my @temp = split(//, $edgeList{$prevCh});
		$shuf_seq .= shift(@temp);
		$edgeList{$prevCh} = join('', @temp);		
	} 
	
	return ($shuf_seq);
}

############################################################################
# Zero-order Markov Model sequence generation
############################################################################
sub genZMMRNASet {

	my ($seq) = @_;
	my $seqLen = length($seq);
	my $shuf_seq;

	#Define the counts of monomers
	my %cnt_monomers = %gl_monomers;
    #Define the frequencies of monomers
	my %freq_monomers = %gl_monomers;
    #Define the cf of monomers
	my %cf_monomers = %gl_monomers;
	
	
    #Compute counts of monomer and dimer 
	for my $i (0..$seqLen-1) {
		my $monomer = substr($seq, $i, 1);
		$cnt_monomers{$monomer}++ if defined $cnt_monomers{$monomer};		
	}	
                 	 
                 	 
	#Compute frequencies of monomer
	foreach my $i (keys(%cnt_monomers)){			
		$freq_monomers{$i} = $cnt_monomers{$i}/$seqLen;				
	}
	
	#compute cf of monomer
	$cf_monomers{'A'} = $freq_monomers{'A'};
	$cf_monomers{'C'} = $cf_monomers{'A'} + $freq_monomers{'C'};
	$cf_monomers{'G'} = $cf_monomers{'C'} + $freq_monomers{'G'};
	$cf_monomers{'U'} = $cf_monomers{'G'} + $freq_monomers{'U'};

	for my $i (0..$seqLen-1) {
		my $randNum = rand;
		foreach my $j (sort keys(%cf_monomers)) {
			if ($randNum <= $cf_monomers{$j}) {
				$shuf_seq .= $j;
				last;
			}
		}
	}		
	
	return($shuf_seq);
}

############################################################################
# First-order Markov Model sequence generation
############################################################################
sub genFMMRNASet {
	
	my ($seq) = @_;
	my $seqLen = length($seq);
	
	#Define the counts of monomers and dimers
	my %cnt_monomers = %gl_monomers;
	my %cnt_dimers = %gl_dimers;
    #Define the frequencies of monomers and dimers
	my %freq_monomers = %gl_monomers;
	my %freq_dimers = %gl_dimers;
    #Define the cf of monomers and dimers
	my %cf_monomers = %gl_monomers;
                   	   	
	my %cumprobMono = %gl_monomers; 
	my %cpdi = %gl_dimers;
	my %cpdi_adj = %gl_dimers;
	my %cupdi = %gl_dimers;

    #Compute counts of monomer and dimer 
	for my $i (0..$seqLen-1) {
		my $monomer = substr($seq, $i, 1);
		$cnt_monomers{$monomer}++ if defined $cnt_monomers{$monomer};
		
		my $dimer = substr($seq, $i, 2);
		$cnt_dimers{$dimer}++ if defined $cnt_dimers{$dimer};	
	}	
	           	 
	#Compute frequencies of monomer
	foreach my $i (keys(%cnt_monomers)) {
		$freq_monomers{$i} = $cnt_monomers{$i}/$seqLen;	
	}
	
	#Compute frequencies of dimers
	foreach my $i (keys(%freq_dimers)) {
		$freq_dimers{$i} = $cnt_dimers{$i}/($seqLen - 1);
	}
		
	#compute cf of monomer
	$cf_monomers{'A'} = $freq_monomers{'A'};
	$cf_monomers{'C'} = $cf_monomers{'A'} + $freq_monomers{'C'};
	$cf_monomers{'G'} = $cf_monomers{'C'} + $freq_monomers{'G'};
	$cf_monomers{'U'} = $cf_monomers{'G'} + $freq_monomers{'U'};

	foreach my $i (sort keys (%cumprobMono)) {
	
		foreach my $j (sort keys (%cumprobMono)){
	    	$cpdi{$i.$j} = $freq_dimers{$i.$j}/$freq_monomers{$i} if ($cnt_monomers{$i} > 0);
	    }    
	    
		foreach my $j (sort keys (%cumprobMono)) {
	    	$cumprobMono{$i} += $cpdi{$i.$j};
	    }

		foreach my $j (sort keys (%cumprobMono)) {
	    	$cpdi_adj{$i.$j} = $cpdi{$i.$j}/$cumprobMono{$i} if ($cumprobMono{$i} > 0);  	
	    }		    

		my $temp = 0;
		foreach my $j (sort keys (%cumprobMono)) {
			$temp += $cpdi_adj{$i.$j};
	    	$cupdi{$i.$j} = $temp;
	    }
	}

	my $randNum = rand;
	my @shuf_seq;
	
	# Select the first random nucleotide
	foreach my $i (sort keys(%cf_monomers)) {	
		if ($randNum <= $cf_monomers{$i}) {
				$shuf_seq[0] = $i;
				last;
		}
	}

	for my $i (1..$seqLen-1) {
		my $randNum = rand;	
		foreach my $j (sort keys(%cf_monomers)) {
			if ($randNum <= $cupdi{$shuf_seq[$i-1].$j}) {
				$shuf_seq[$i] = $j;
				last;
			}
		}
	}
	
	return(join("", @shuf_seq));
}

exit;

