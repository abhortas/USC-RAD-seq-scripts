#! /usr/bin/perl

# This script is make_forward.pl

# Author: Andrés Blanco Hortas (andres.blanco.hortas@usc.es)

####################################################################################################################################################################################
# This script is intended to flip sequences with palindromic recognition site (e.g. AlfI), to have them all orientated in the same strand.
# Use this script when there is no reference genome.

# Fastq format is required as input.
# This script require the following free tools: "Bowtie 1" (https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.0/)
# and "CD-HIT-EST" (http://weizhongli-lab.org/cd-hit/).

# we recommend checking the parameters. 

# References:


# BOWTIE 1:

# Langmead, B., Trapnell, C., Pop, M. et al. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol 10, R25 (2009). https://doi.org/10.1186/gb-2009-10-3-r25
# Langmead, B. (2010), Aligning Short Sequencing Reads with Bowtie. Curr. Protoc. Bioinform., 32: 11.7.1-11.7.14. doi:10.1002/0471250953.bi1107s32

# CD-HIT:

# Fu L, Niu B, Zhu Z, Wu S, Li W. CD-HIT: Accelerated for clustering the next-generation sequencing data. Bioinformatics. 2012;28:3150–2. https://doi.org/10.1093/bioinformatics/bts565.
# Li W, Godzik A. Cd-hit: A fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics. 2006;22:1658–9.https://doi.org/10.1093/bioinformatics/btl158.


##########################################################################################

use strict;
use warnings;
use Getopt::Long;



my ($dir,$catalog,$c,$cdhit,$output,$bowtie,$reverse,$g,$v);
my ($line,$archive,$i,$elem,$element,$help,$k);
my @archive;
$c=0.916;
$g=1;
$v=3;

################################### Get options ##########################################
GetOptions(	"dir|d=s"  		=> \$dir,
			"catalog|ct=s" => \$catalog,
			"c|c=s"			=> \$c,
			"g|g=s"			=> \$g,
			"v|v=s"			=> \$v,
			"output|o=s"	=> \$output,
			"help|h"     	=> \$help,) or die "Try 'make_forward.pl --help' for more information\n"; 	
##########################################################################################

################################## Check the options #####################################
if($help){
	print "Use: perl make_forward.pl [OPTION]...\n";
	print "Script to flip sequences with palindromic recognition site in the same strand\n";
	print "\n";
	print "Options:\n";
	print "\t-d,  --dir\t\tPath to the directory with input Fastq\n";
	print "\t-ct, --catalog\t\tarchives catalog (directory or archive with the names).The same directory as -d\n";
	print "\t-c,  --c\t\tOption -c from ch-hit-est (i.e. sequence identity threshold).\n";
	print "\t-g,  --g\t\tOption -g from ch-hit-est (best match; g=1 default).\n";
	print "\t-v,  --v\t\tOption -v from Bowtie 1 (report alignments with at most <int> mismatches; v=3 default).\n";
	print "\t-o,  --output\t\tPath to the directory for output files.\n";
	print "\t-h,  --help\t\tShow this help message and exit.\n";
	print "\n";
	exit;
}

unless($dir){
	print "Error: The option --dir is necessary.\n";
	exit;
}

unless(-d $dir){
	print "Error: The specified directory with input files do not exist.\n";
	exit;
}

unless($c){
	print "Error: The option --catalog is necessary.\n";
	exit;
}
#######################################################


################################# subroutine fastq to fasta
sub fast{

	my ($in,$out)=@_;
	my ($k,$line);
	
	$k=0;
	open(FILE,"$in");
	open(TEXT,">$out");
	while($line = <FILE>){
		if($k == 1){
			print TEXT $line;
			$k=0;
		}
		if($line =~ /^@/){
			$k=1;
			$line =~ s/@/>/;
			print TEXT $line;
		}
	}
	close(FILE);
	close(TEXT);
}
#################################

################################# subroutine reverse
sub reverse_fq{

	my ($data1,$out)=@_;
	my ($line,$i,$j,$k,$new_seq,$new_cal);
	my @a;

	open(FILE1,"$data1");
	open(TEXT,">$out");
	while($line = <FILE1>){
		unless($line =~ /^@/){
			chomp($line);
			@a=split(/\t/,$line);
			print TEXT "@","$a[0]\n$a[9]\n+\n$a[10]\n";
		}
	}
	close(FILE1);
	close(TEXT);
}

################################# A directory with the catalog files is created
mkdir("$output/catalog");
mkdir("$output/references");
if(-f $catalog){
	open(FILE,"$catalog");
	$i=-1;
	while($line = <FILE>){
		$i++;
		$line =~ s/\r\n/\n/;
		chomp($line);
		$archive[$i]=$line;
	}
	close(FILE);
	chdir $dir;
	foreach $elem(<*>){
		$k=0;
		for($i=0; $i<=$#archive; $i++){
			if($elem eq $archive[$i]){
				$k=1;
				last;
			}
		}
		if($k == 1){
			$element=$elem;
			$element =~ s/.fastq//;
			& fast ($elem,"$output/catalog/$element\.fasta");
			system("cd-hit-est -i $output/catalog/$element\.fasta -c 1 -g $g -T 0 -M 0 -d 0 -o $output/references/$element\.fasta");
		}
	}
}

elsif(-d $catalog){
	chdir $catalog;
	foreach $elem(<*>){
		$element=$elem;
		$element =~ s/.fastq//;
		& fast ($elem,"$output/catalog/$element\.fasta");
		system("cd-hit-est -i $output/catalog/$element\.fasta -c 1 -g $g -T 0 -M 0 -d 0 -o $output/references/$element\.fasta");		
	}
}

mkdir("$output/reference");
system("cat $output/references/*.fasta > $output/reference/reference.fasta");
system("cd-hit-est -i $output/reference/reference.fasta -c $c -g $g -T 0 -M 0 -d 0 -o $output/reference.fasta");
mkdir("$output/index");
system("bowtie-build $output/reference.fasta $output/index/reference");
###################################

################################### Align each file
chdir $dir;
mkdir("$output/aligned");
mkdir("$output/definitives");
foreach $elem(<*>){
	$element=$elem;
	$element =~ s/.fastq//;
	system("bowtie --best --strata -m 2 -k 1 -p 24 --sam -q -v $v --chunkmbs 50000 $output/index/reference $elem > $output/aligned/$element.sam");
	& reverse_fq ("$output/aligned/$element.sam","$output/definitives/$elem");
}
#####################################	



