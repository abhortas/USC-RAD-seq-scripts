#! usr/bin/perl

# This script is new_sam.pl

# Author: Andrés Blanco Hortas (andres.blanco.hortas@usc.es)

####################################################################################################################################################################################

# This script is for reference-genome approach, for single-end sequencing.
# To modify input sam files with aligned reads. To convert all aligned reverse reads (i.e. SAM flag "16") into aligned forward reads (i.e. SAM flag "0").
# Minor modifications would be necessary with paired-end sequencing.

# To decoding SAM flags you can use: https://broadinstitute.github.io/picard/explain-flags.html

# Broad Insitute. Picard Tools. http://broadinstitute.github.io/picard/. Accessed 1 July 2020.

####################################################################################################################################################################################

use strict;
use warnings;
use Getopt::Long;

my ($input,$out,$line,$resto,$elem,$help);
my @a;

#################################### Get options #########################################
GetOptions(	"input|i=s" 	 => \$input,
			"output|o=s"	 => \$out,
			"help|h"  	     => \$help,) or die "Try 'perl new_sam.pl --help' for more information\n"; 	
##########################################################################################

#################################### Check the arguments #################################
if($help){
	print "Use: perl new_sam.pl [OPTION]...\n";
	print "Script to create a new .sam file with every alignment in forward strand (i.e. with SAM flag 0)\n";
	print "\n";
	print "Options:\n";
	print "\t-i, --input\t\tAbsolute path to the directory with the input files (.sam).\n";
	print "\t-o, --output\t\tAbsolute path to the difectory for output files.\n";
	print "\t-h, --help\t\tShow this help message and exit.\n";
	print "\n";
	exit;
}
unless($input){
	print "Error: The option --input is necessary.\n";
	exit;
}

unless($out){
	print "Error: The option --output is necessary.\n";
	exit;
}

unless(-d $input){
	print "Error: The specified directory for the input files does not exist.\n";
	exit;
}
unless(-d $out){
	print "Error: The specified directory for the output files does not exist.\n";
	exit;
}
##########################################################################################

###################################### Program ##########################################
chdir $input;
foreach $elem(<*>){
	open(FILE,"$elem");
	open(TEXT,">$out/$elem");
	while($line = <FILE>){
		if($line =~ /^@/){
			print TEXT $line;
		}
		else{
			chomp($line);
			unless($line =~ /^$/){
				@a=split(/\t/,$line);
				if($a[1] == 16){
					$resto=join("\t",@a[2..$#a]);
					print TEXT "$a[0]\t0\t$resto\n";				
				}
				else{
					print TEXT "$line\n";
				}
			}
		}
	}
	close(FILE);
	close(TEXT);
	print "File $elem modified.\n";
}
exit;
##########################################################################################
			
