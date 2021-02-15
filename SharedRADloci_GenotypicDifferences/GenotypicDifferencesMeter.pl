#! /usr/bin/perl

#This script is GenotypicDifferencesMeter.pl

#Author: Andrés Blanco Hortas (andres.blanco.hortas@usc.es)

##########################################################################################

# References:

#About Genepop format: https://genepop.curtin.edu.au/help_input.html

##########################################################################################


use strict;
use warnings;
use Getopt::Long;

my ($shared,$genepop1,$genepop2,$line,$out,$help,$g11,$g12,$g21,$g22,$k,$i,$j);
my ($homo_homo,$het_mis,$mis_het,$homo_mis,$mis_homo,$het_homo,$homo_het,$equal);
my (@a,@ind,@loci);
my (%geno1,%geno2);


################################### Get options ##########################################
GetOptions(	"shared|s=s"  		=> \$shared,
			"genepop1|g1=s" 	=> \$genepop1,
			"genepop2|g2=s"		=> \$genepop2,
			"output|o=s"		=> \$out,
			"help|h"     		=> \$help,) or die "Try 'GenotypicDifferencesMeter.pl --help' for more information about options\n";

##########################################################################################

################################## Check the options #####################################
if($help){
	print "Use: perl GenotypicDifferencesMeter.pl [OPTIONS]...\n\n";
	print "\n";
	print "This script use 3 input files:\n";
	print "A shared SNPs index (output from SharedExclusiveSNPs.pl script) with three columns: SNP pair 1<tab>SNP pair 2<tab>strand (i.e. +/-).\n";
	print "Two genepops from different building-loci pipelines. The genepop format is the same as the STACKS output:\n"; 	
	print "\n";
	print "Genepop format: separate loci with commas on the same line.\n";
	print "Line 2 : the name of the Nth locus.\n";
	print "Line 3 : Type the word Pop, or pop.\n";
	print "Line 4 : And example is given below:\n";
	print "ind#001fem, 0101 0202 0000 0104\n";
	print "\n";
	print "Options:\n";
	print "\t-s,   --shared\t\tShared SNPs file (output from SharedExclusiveSNPs.pl script).\n";
	print "\t-g1,  --genepop1\tFirst genepop file (file corresponding to the SNPs of first column of Shared SNPs index).\n";
	print "\t-g2,  --genepop2\tSecond genepop file (file corresponding to the SNPs of second column of Shared SNPs index).\n";
	print "\t-o,   --output\t\tOutput file path\n";
	print "\t-h,   --help\t\tShow this help message and exit.\n";
	print "\n";
	exit;
}

unless($shared){
	print "Error: The option --shared is necessary.\n";
	exit;
}

unless(-e $shared){
	print "Error: The file with the Shared SNPs does not exist.\n";
	exit;
}

unless($genepop1){
	print "Error: The option --genepop1 is necessary.\n";
	exit;
}

unless(-e $genepop1){
	print "Error: The genepop1 file does not exist.\n";
	exit;
}

unless($genepop2){
	print "Error: The option --genepop2 is necessary.\n";
	exit;
}

unless(-e $genepop2){
	print "Error: The genepop2 file does not exist.\n";
	exit;
}

unless($out){
	print "Error: The option --output is necessary.\n";
	exit;
}
##########################################################################################

$homo_homo=$het_mis=$mis_het=$homo_mis=$mis_homo=$het_homo=$homo_het=$equal=0;


### Read genotypes from first genepop file
$k=0;
$j=-1;
open(FILE,"$genepop1");
while($line = <FILE>){
	$k++;
	if($k == 2){
		chomp($line);
		@loci=split(/,/,$line);
	}
	elsif($k > 2){
		unless(($line =~ /^pop/) | ($line =~ /^Pop/)){
			chomp($line);
			@a=split(/\t/,$line);
			$j++;
			$a[0] =~ s/,//;
			$ind[$j]=$a[0];
			for($i=1; $i<=$#a; $i++){
				$geno1{$loci[$i-1]}{$a[0]}=$a[$i];
			}
		}
	}
}
close(FILE);


### Read genotypes from second genepop file
$k=0;
$j=-1;
open(FILE,"$genepop2");
while($line = <FILE>){
	$k++;
	if($k == 2){
		chomp($line);
		@loci=split(/,/,$line);
	}
	elsif($k > 2){
		unless(($line =~ /^pop/) | ($line =~ /^Pop/)){
			chomp($line);
			@a=split(/\t/,$line);
			$j++;
			$a[0] =~ s/,//;
			for($i=1; $i<=$#a; $i++){
				$geno2{$loci[$i-1]}{$a[0]}=$a[$i];
			}
		}
	}
}
close(FILE);


# Calculate differences
open(FILE,"$shared");
while($line = <FILE>){
	chomp($line);
	@a=split(/\t/,$line);
	for($i=0; $i<=$#ind; $i++){
		if($a[2] eq "-"){
			$geno2{$a[1]}{$ind[$i]} =~ tr/1234/4321/;
		}
		$g11=substr($geno1{$a[0]}{$ind[$i]},0,2);
		$g12=substr($geno1{$a[0]}{$ind[$i]},2,2);
		$g21=substr($geno2{$a[1]}{$ind[$i]},0,2);
		$g22=substr($geno2{$a[1]}{$ind[$i]},2,2);
		if($geno1{$a[0]}{$ind[$i]} eq "0000"){
			if($geno2{$a[1]}{$ind[$i]} eq "0000"){
				$equal++;
			}
			else{
				if($g21 eq $g22){
					$mis_homo++;
				}
				else{
					$mis_het++;
				}
			}
		}
		else{
			if($g11 eq $g12){
				if($geno2{$a[1]}{$ind[$i]} eq "0000"){
					$homo_mis++;
				}
				else{
					if($g21 eq $g22){
						if($g11 eq $g22){
							$equal++;
						}
						else{
							$homo_homo++;
						}
					}
					else{
						$homo_het++;
					}
				}
			}
			else{
				if($geno2{$a[1]}{$ind[$i]} eq "0000"){
					$het_mis++;
				}					
				else{
					if($g21 eq $g22){
						$het_homo++;
					}
					else{
						$equal++;
					}
				}
			}
		}
	}
}
		

# write results
open(TEXT,">$out");
print TEXT "Equal genotypes:\t$equal\n";
print TEXT "Homo1-Homo2:\t$homo_homo\n";
print TEXT "Homo-Het:\t$homo_het\n";
print TEXT "Homo-Mis:\t$homo_mis\n";
print TEXT "Het-Homo:\t$het_homo\n";
print TEXT "Het-Mis:\t$het_mis\n";
print TEXT "Mis-Homo:\t$mis_homo\n";
print TEXT "Mis-Het:\t$mis_het\n";
				
