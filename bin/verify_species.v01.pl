#! /usr/bin/perl

# Program to compare the most abundant species detected in the illumina and nanopore reads with the LIMS record for rthe specified sample
# Information from the reads is in one-line format species1:abundance1;....;speciesN:abundanceN;  
# see the SAMPLE.[illumina/ont].bracken_abundance.tophits  files

# LIMS does not work from computational nodes. Requires internet connection, only available from login nodes

# Fernando Cruz, 2019-07-12

use strict;
use Getopt::Long;
use File::Basename;

# define variables

my $illumina=""; 
my $nanopore=""; 
my $seq_sample="";

my @fields=();
my $lims_species="";

my $illumina_line="";
my @bracken_illumina=();
my $illumina_species="";

my $nanopore_line="";
my @bracken_nanopore=();
my $nanopore_species="";

my $output="";
my $warning="";

my $agree="False";
my $verify="False";

# get these two files with the bracken results
GetOptions(
           'i:s'  => \$illumina,
           'n:s'  => \$nanopore,
           'l:s'  => \$lims_species
);


if ( ($illumina eq "") || ($nanopore eq "") || ($lims_species eq "") ) {
       die print "provide bracken abundance for illumina and nanopore/ONT reads! Also requires the value of lims species as parameter \n\tUSAGE: verify_species.v01.pl \-i illumina.bracken_abundance.tophits \-n ONT.bracken_abundance.tophits \-l \"lims_species\"\n\t\t\tNote that double tilde matters...\n"; 
}

my ($illumina_name,$Ipath,$Iext) = fileparse($illumina,qw(\.illumina\.bracken_abundance\.tophits));
my ($nanopore_name,$Npath,$Next) = fileparse($nanopore,qw(\.ont\.bracken_abundance\.tophits));

if ($illumina_name eq $nanopore_name){
    $seq_sample=$illumina_name;
   # print "sample name is $sample\n";
}
else{
    die print "sample name for nanopore and illumina reads differ!\nPlease check them\:\n$nanopore_name\t$illumina_name\n";
}

# parsing illumina tophit
my $count=0;
open (IN, "$illumina") or die print "cannot open $illumina file\n";
while (<IN>){
    next unless ($count==0);#get just first line
    $illumina_line=$_;
    chomp; 
    #print "$_\n";
    @bracken_illumina=split/\t/,$_;
    $illumina_species=$bracken_illumina[1];
    @bracken_illumina=split/\:/,$illumina_species;
    $illumina_species=$bracken_illumina[0];
    #print "illumina $illumina_species\n";
    $count++;
}
close (IN) or die print "cannot close $illumina file\n";


#parsing nanopore tophit
$count=0;
open (IN2, "$nanopore") or die print "cannot open $nanopore file\n";

while (<IN2>){
  
   next unless ($count==0);#get just first line
    $nanopore_line=$_;
   chomp; 
   #print "$_\n";
   @bracken_nanopore=split/\t/,$_;
   $nanopore_species=$bracken_nanopore[1];
   @bracken_nanopore=split/\:/,$nanopore_species;
   $nanopore_species=$bracken_nanopore[0];
   #print "nanopore $nanopore_species\n";
   $count++;
}
close (IN2) or die print "cannot close $nanopore file\n";

#print "species name illumina\: $illumina_species nanopore\: $nanopore_species LIMS\: $lims_species\n";


if ("$illumina_species" eq "$nanopore_species") { 
    $agree="True"; 
}

if ( ("$illumina_species" eq "$lims_species") && ("$nanopore_species" eq "$lims_species") ) { 
    $verify="True";
}

if ($lims_species eq ""){
   $verify="NA";
}

# send output to a file
$output="$seq_sample\."."species_verification\.out";


open (OUT, ">$output") or die "cannot open $output in current directory\n";

print OUT "Illumina_sp\tONT_sp\tLIMS_sp\tSeq_Agreement\tVerified\n";
print OUT "$illumina_species\t$nanopore_species\t$lims_species\t$agree\t$verify\n";

close (OUT) or die "cannot close $output in current directory\n";



# Write to standard output this message when species don't match. This will go into the log
if ( ($agree eq "False") || ($verify eq "False") ) {     
    print "\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\SPECIES MISMATCH\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\n"; 
    print "There is some disagreement between species name in bracken results and/or the LIMS.\n";
    print "Please check below\:\n";
    print "$seq_sample Illumina Bracken abundances\:\n\t$illumina_line\n";
    print "$seq_sample ONT Bracken abundances\:\n\t$nanopore_line\n";
    print "$seq_sample one LIMS record\:\n\t$lims_species\n";
    print "\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\n";
}
else {
    print "$seq_sample sequencing data agree with the LIMS\n";
}
