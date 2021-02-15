#! /usr/bin/perl

#This script is SharedExclusiveSNPs.pl

#Author: Andrés Blanco Hortas (andres.blanco.hortas@usc.es)

##########################################################################################

# References for the mentioned software in this perl script

# Building-loci pipelines:

# STA: Catchen J, Hohenlohe PA, Bassham S, Amores A, Cresko WA. Stacks: An analysis tool set for population genomics. Mol Ecol. 2013;22:3124–40. https://doi.org/10.1111/mec.12354.
# STA: Catchen JM, Amores A, Hohenlohe P, Cresko W, Postlethwait JH. Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. G3. 2011;1:171–82.https://doi.org/10.1534/g3.111.000240.
# ALT: Wang S, Meyer E, McKay JK, Matz M V. 2b-RAD: a simple and flexible method for genome-wide genotyping. Nat Methods. 2012;9:808–10. https://doi.org/10.1038/nmeth.2023.

# CD-HIT:

# Fu L, Niu B, Zhu Z, Wu S, Li W. CD-HIT: Accelerated for clustering the next-generation sequencing data. Bioinformatics. 2012;28:3150–2. https://doi.org/10.1093/bioinformatics/bts565.
# Li W, Godzik A. Cd-hit: A fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics. 2006;22:1658–9.https://doi.org/10.1093/bioinformatics/btl158.

##########################################################################################

use strict;
use warnings;
use Getopt::Long;


my ($data,$meyer,$stacks,$out,$line,$k,$rad1,$rad2,$or1,$or2,$length,$ind,$elem,$or,$count,$pos,$i,$j,$help,$script,$ultimo,$primer,$p);
my ($mexcl,$sexcl,$comp);
my (@a,@b,@ps,@pm,@rads);
my (%meyer1,%meyer2,%meyer3,%stacks1,%stacks2,%stacks3);


################################### Get options ##########################################
GetOptions(	"stacks|s=s"  	=> \$stacks,
			"meyer|m=s" 	=> \$meyer,
			"cluster|c=s"	=> \$data,
			"output|o=s"	=> \$out,
			"length|l=s"		=> \$length,
			"help|h"     	=> \$help,) or die "Try 'perl shared_exclusive_SNPs_v2.pl --help' for more information about options\n";

##########################################################################################

################################## Check the options #####################################
if($help){
	print "Use: perl SharedExclusiveSNPs.pl [OPTIONS]...\n";
	print "It is necessary that the first sequences of the fasta input for CD-HIT belong to STACKS.\n";
	print "\n";
	print "This script use the .clstr output from CD-HIT. Furthermore, it uses two more input files consistent in SNPs index from each pipeline (i.e. -s and -m options).
		This script was designed to obtain shared and exclusive SNPs from STACKS 2 and Meyer's 2b-RAD v2.1 pipeline.
		For comparisons between different bulding-loci pipelines, minor modifications are necessary.\n";
	print "\n";
	print "Only would be considered shared SNPs those from CD-HIT clusters composed by only one RAD-loci per pipeline.
		It would be considered exclusive SNPs those from CH-HIT clusters composed by RAD-loci belonging to the same building-loci pipeline.
		Uniquely SNPs which RAD-loci cluster together and share position will be considered as shared SNPs.
		You would have to check if both pipelines under comparison have the same counting system (i.e. if the first position is 0 or 1).\n";
	print "\n";
	print "Options:\n";
	print "\t-s,  --stacks\t\tFile with the STACKS SNP index.\n";
	print "\t-m,  --meyer\t\tFile with the Meyer SNP index.\n";
	print "\t-c,  --cluster\t\tCluster file from CD-HIT.\n";
	print "\t-o,  --output\t\toutput path to write result files\n";
	print "\t-l,  --long\t\tLength of the compared sequences.\n";
	print "\t-h,  --help\t\tShow this help message and exit.\n";
	print "\n";
	exit;
}

unless($stacks){
	print "Error: The option --stacks is necessary.\n";
	exit;
}

unless(-e $stacks){
	print "Error: The file with the STACKS SNP index does not exist.\n";
	exit;
}

unless($meyer){
	print "Error: The option --meyer is necessary.\n";
	exit;
}

unless(-e $meyer){
	print "Error: The file with the Meyer SNP index does not exist.\n";
	exit;
}

unless($data){
	print "Error: The option --cluster is necessary.\n";
	exit;
}

unless(-e $data){
	print "Error: The clusters file from CD-HIT does not exist.\n";
	exit;
}

unless($out){
	print "Error: The option --output is necessary.\n";
	exit;
}

unless(-e $out){
	print "Error: The specified directory for the output files does not exist.\n";
	exit;
}

unless($length){
	print "Error: The option --long is necessary.\n";
	exit;
}
##########################################################################################



##########################################################################################



open(FILE,"$meyer");
while($line = <FILE>){
	chomp($line);
	@a=split(/_/,$line);
	if(exists $meyer1{$a[0]}){
		if(exists $meyer2{$a[0]}){
			if(exists $meyer3{$a[0]}){
				print "ERROR: There are more than 3 SNPs at least into one meyer rad, and this is not allowed.\n";	
				exit;
			}
			else{
				$meyer3{$a[0]}=$a[1];
			}
		}
		else{
			$meyer2{$a[0]}=$a[1];
		}
	}
	else{
		$meyer1{$a[0]}=$a[1];
	}
}
close(FILE);


$k=-1;
open(FILE,"$stacks");
while($line = <FILE>){
	chomp($line);
	@a=split(/_/,$line);
	if(exists $stacks1{$a[0]}){
		if(exists $stacks2{$a[0]}){
			if(exists $stacks3{$a[0]}){
				print "ERROR: There are more than 3 SNPs at least into one stacks rad, and this is not allowed.\n";	
				exit;			
			}
			else{
				$stacks3{$a[0]}=$a[1];
			}
		}
		else{
			$stacks2{$a[0]}=$a[1];
		}
	}
	else{
		$stacks1{$a[0]}=$a[1];
	}
}
close(FILE);

$ultimo="";
$comp=0;
$sexcl=0;
$mexcl=0;
$k=0;
open(TEXT,">$out/shared_SNPs.txt");
open(UNIS,">$out/exclusive_stacks.txt");
open(UNIM,">$out/exclusive_meyer.txt");
open(FILE,"$data");
while($line = <FILE>){
	if($line =~ /^>/){
		$k++;
		if($k > 1){
			if($ultimo eq ""){
				if($primer eq "meyer"){
					print UNIM "$rads[0]\_",$meyer1{$rads[0]},"\n";
					$mexcl++;
				}
				else{
					print UNIS "$rads[0]\_",$stacks1{$rads[0]},"\n";
					$sexcl++;
				}
			}
			else{
				if($primer eq $ultimo){
					if($primer eq "meyer"){
						foreach $elem(@rads){
							print UNIM "$elem\_",$meyer1{$elem},"\n";
							$mexcl++;
						}
					}
					else{
						foreach $elem(@rads){
							print UNIS "$elem\_",$stacks1{$elem},"\n";
							$sexcl++;
						}
					}
				}
			}
			undef(@rads);										
			$ultimo="";
			$primer="";
			if(($count == 1) & ($i == 2)){
				@ps=("x","x","x");
				@pm=("y","y","y");
				$ps[0]=$stacks1{$rad1};
				if(exists $stacks2{$rad1}){
					$ps[1]=$stacks2{$rad1};
				}
				if(exists $stacks3{$rad1}){
					$ps[2]=$stacks3{$rad1};
				}
				if($or eq "+"){
					$pm[0]=$meyer1{$rad2};
					if(exists $meyer2{$rad2}){
						$pm[1]=$meyer2{$rad2};
					}
					if(exists $meyer3{$rad2}){
						$pm[2]=$meyer3{$rad2};
					}
				}
				else{
					$pm[0]=$length-$meyer1{$rad2}+1;
					if(exists $meyer2{$rad2}){
						$pm[1]=$length-$meyer2{$rad2}+1;
					}
					if(exists $meyer3{$rad2}){
						$pm[2]=$length-$meyer3{$rad2}+1;
					}					
				}
				$j=0;
				foreach $elem(@ps){
					foreach $ind(@pm){
						if($elem eq $ind){
							if($or eq "+"){
								print TEXT "$rad1\_$elem\t$rad2\_$ind\t+\n";
							}
							else{
								$pos=$length-$ind+1;
								print TEXT "$rad1\_$elem\t$rad2\_$pos\t-\n";
							}
							$j=1;
							$comp++;
						}
					}
					if($j == 1){
						last;
					}
				}
			}
		}		
		$i=0;	
		$p=0;		
	}
	else{
		@a=split(/\t/,$line);
		$count=$a[0];
		if($count == 0){
			if($a[1] =~ /denovoLocus/){
				$a[1] =~ />(.+?)\./;
				$rads[$p]=$1;
				$primer="meyer";
				$p++;
			}
			else{
				$a[1] =~ />(.+?)\./;
				$rad1=$1;
				$rads[$p]=$rad1;
				$primer="stacks";
				$p++;
				$i++;
			}
		}
		elsif($count == 1){
			if($a[1] =~ /denovoLocus/){
				$a[1] =~ />(.+?)\./;
				$rad2=$1;
				$rads[$p]=$rad2;	
				$a[1] =~ /at (.+?)\//;	
				$or=$1;	
				$ultimo="meyer";
				$i++;
				$p++;
			}
			else{
				$a[1] =~ />(.+?)\./;
				$rads[$p]=$1;
				$ultimo="stacks";
				$p++;
			}
		}
		else{
			if($a[1] =~ /denovoLocus/){
				$a[1] =~ />(.+?)\./;
				$rad2=$1;
				$rads[$p]=$rad2;	
				$ultimo="meyer";
				$p++;
			}
			else{
				$a[1] =~ />(.+?)\./;
				$rads[$p]=$1;
				$ultimo="stacks";
				$p++;
			}
		}				
	}
}


# Ultimo exclusivo
if($ultimo eq ""){
	if($primer eq "meyer"){
		print UNIM "$rads[0]\_",$meyer1{$rads[0]},"\n";
		$mexcl++;
	}
	else{
		print UNIS "$rads[0]\_",$stacks1{$rads[0]},"\n";
		$sexcl++;
	}
}
else{
	if($primer eq $ultimo){
		if($primer eq "meyer"){
			foreach $elem(@rads){
				print UNIM "$elem\_",$meyer1{$elem},"\n";
				$mexcl++;
			}
		}
		else{
			foreach $elem(@rads){
				print UNIS "$elem\_",$stacks1{$elem},"\n";
				$sexcl++;
			}
		}
	}
}


# Ultimo compartido
if(($count == 1) & ($i == 2)){
	@ps=("x","x","x");
	@pm=("y","y","y");
	$ps[0]=$stacks1{$rad1};
	if(exists $stacks2{$rad1}){
		$ps[1]=$stacks2{$rad1};
	}
	if(exists $stacks3{$rad1}){
		$ps[2]=$stacks3{$rad1};
	}
	if($or eq "+"){
		$pm[0]=$meyer1{$rad2};
		if(exists $meyer2{$rad2}){
			$pm[1]=$meyer2{$rad2};
		}
		if(exists $meyer3{$rad2}){
			$pm[2]=$meyer3{$rad2};
		}
	}
	else{
		$pm[0]=$length-$meyer1{$rad2}+1;
		if(exists $meyer2{$rad2}){
			$pm[1]=$length-$meyer2{$rad2}+1;
		}
		if(exists $meyer3{$rad2}){
			$pm[2]=$length-$meyer3{$rad2}+1;
		}					
	}
	$j=0;
	foreach $elem(@ps){
		foreach $ind(@pm){
			if($elem eq $ind){
				if($or eq "+"){
					print TEXT "$rad1\_$elem\t$rad2\_$ind\t+\n";
				}
				else{
					$pos=$length-$ind+1;
					print TEXT "$rad1\_$elem\t$rad2\_$pos\t-\n";
				}
				$j=1;
				$comp++;
			}
		}
		if($j == 1){
			last;
		}
	}
}

print "Shared SNPs between pipelines: $comp\n";
print "Exclusive SNPs from STACKS: $sexcl\n";
print "Exclusive SNPs from Meyer: $mexcl\n";


				 



