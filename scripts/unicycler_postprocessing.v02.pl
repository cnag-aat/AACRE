#! /usr/bin/env perl

# Takes Unicycler format and produces a table (CRE project)
#       includes the depth (inferred Stoichiometry) as a decimal numer
#       prints out assembly with new sequenceID

# Fernando Cruz, 2019-10-11
use strict;
use Getopt::Long;
use File::Basename;

# define variables

my $sample='';
my $file=''; 
my $assembler='';
my $pipeversion='';

my @line;# array for parsing elements
my $seqid;
my $contig;
my $length;
my $depth;
my $circular;

# new variables
my $assembly_name="none";
my $total_scaffolds=0;
my $total_circular=0;
my $circularity_ratio=0;
my $assembly_length=0;
my $max_scaffold_length=0;
my $scaffolds_2kb_or_shorter=0;

GetOptions(
           'f|file:s'  => \$file,
           's|sample:s' => \$sample,
           'a|assembler:s' => \$assembler,
           'p|pipeversion:s'=> \$pipeversion
);

#require all options
if ( ($file eq "") || ($sample eq "") || ($assembler eq "") || ($pipeversion eq "") ){
       die print "mandatory to provide unicycler assembly filename, sample barcode, assembler used and snakemake pipeline version\. \nUSAGE\: unicycler2table.pl \-f filename \-s sample_barcode \-a assembler_used \-p pipeline_version\n"; 
}

open (TABLE_ASM, ">tables\/$sample\.$pipeversion\.assembly\.tbl") or die print "cannot open table file \n";
open (TABLE_SCF, ">tables\/$sample\.$pipeversion\.scaffolds\.tbl") or die print "cannot open table file \n";
open (NEWFASTA, ">$sample\.$pipeversion\.assembly.fasta") or die print "cannot open new fasta file \n";

#AH0315.v1.assembly.fasta
  #print table header
   print TABLE_SCF "assembly\tscaffold\tscaffold_length\tdepth\tcircular\n"; 
   print TABLE_ASM "sample_barcode\tassembly\ttotal_scaffolds\tcircular_scaffolds\tcircularity_ratio\tscaffolds_2kb_or_shorter\tassembly_length\tmax_scaffold_length\tassembler\n";
 
# reading FASTA
   open (IN, "$file") or die "print cannot open $file\n";

   

   while(<IN>)
   {

      #parse FASTA file descriptions/headers
      chomp;
      if ($_ =~ /^>/) { # print table and modify sequence header
          @line=split /\s+/, $_;

          $contig=$line[0];
          $contig =~ s/\>//g;
          $total_scaffolds+=1;
 
          $assembly_name=$sample."\_".$pipeversion; # pipeversion should already include the "v" 
          $seqid=$sample."\_".$pipeversion."\_c".$contig; # pipeversion should already include the "v" 
 
          $length=$line[1];   
          $length =~ s/length\=//g;# replace to get rid of length
	  if ($contig == 1){
	      $max_scaffold_length=$length;
	  }
          if($length <= 2000){
	      $scaffolds_2kb_or_shorter+=1;
	  }
	  $assembly_length+=$length;# cumulative length of the assembly

          $depth=$line[2];   
          $depth =~ s/depth\=//g;# replace to get rid of depth and equal
          $depth =~ s/x//g;# replace to get rid of the "x" to get just a decimal number

             if (exists $line[3]) {
                $circular=$line[3]; 
                $circular =~ s/circular\=//g;# replace to get rid of length
                $circular =~ s/true/True/g;# replace to put in uppercase
             }
 
             else {
                $circular='False'
	     }
  
        if ($circular eq "True") {
             $total_circular += 1;
        }       
          print TABLE_SCF "$assembly_name\t$seqid\t$length\t$depth\t$circular\n";#scaffolds
          print NEWFASTA "\>$seqid\n";

      }# print new fasta header & table entry
    
      else {
       #print sequence
	  print NEWFASTA "$_\n";
      }
     
   }# END OF INPUT FASTA FILE 

#COMPUTE RATIO OF CONNECTED COMPONENTS (circular scaffolds-or-contigs)
$circularity_ratio= $total_circular/$total_scaffolds;
print TABLE_ASM "$sample\t$assembly_name\t$total_scaffolds\t$total_circular\t$circularity_ratio\t$scaffolds_2kb_or_shorter\t$assembly_length\t$max_scaffold_length\t$assembler\n";

close (IN) or die "print cannot close $file\n";

close (NEWFASTA) or die "print cannot close $sample\.$pipeversion\.assembly\.fasta\n";# larger files closed earlier

close (TABLE_SCF) or die "print cannot close $sample\.$pipeversion\.scaffolds\.tbl\n";

close (TABLE_ASM) or die "print cannot close $sample\.$pipeversion\.assembly\.tbl\n";
