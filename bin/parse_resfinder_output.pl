#!/usr/bin/perl

#Author: Jessica Gomez-Garrido, CNAG-CRG.
#Contact email: jessica.gomez@cnag.crg.eu
use strict;
use warnings;

my $dir = $ARGV[0];
my @ab = ("trimethoprim","tetracycline","sulphonamide","rifampicin","quinolone","phenicol","oxazolidinone","nitroimidazole","macrolide","glycopeptide","fusidicacid", "fosfomycin","colistin","aminoglycoside","beta-lactam");
my $len = scalar @ab; 
my $i = 0;
while ($i < $len){
  my $file = $dir . "/" . "out_$ab[$i].xml";
  open IN, "<", "$file";
  my $id;
  my %blast;
  my $hit;
  while (<IN>){
    chomp;
    my $target;
    if ($_ =~ m/Iteration_query-def\>([^\<]+)/){
      my @id = split / /, $1;
      $id = $id[0];
      #print "$id\n";
    }
    elsif ($_ =~ m/<Hit_num>([^\<]+)/){
      $hit = $1;
      #print "$hit\n";
    }
    elsif ($_ =~ m/Hit_accession\>([^\<]+)/){
      $target= $1;
      $blast{$id}{$hit}->{target}=$target;
    }
    elsif ($_ =~ m/\<Hit_def\>([^\<]+)/){
      my $details = $1;
      $blast{$id}{$hit}->{tstart}=-1;
      $blast{$id}{$hit}->{tend}=-1;
      $blast{$id}{$hit}->{definition} = $details;
     # print "$sp\t$def\t$tg\n";
    }
    elsif ($_ =~ m/<Hsp_evalue>([^\<]+)/){
      $blast{$id}{$hit}->{evalue} = $1;
    }
    elsif ($_ =~m/<Hsp_hit-from>([^\<]+)/){
      $blast{$id}{$hit}->{tstart} = $1 if ($1 < $blast{$id}{$hit}->{tstart} || $blast{$id}{$hit}->{tstart} < 0);
    }
    elsif ($_ =~m/<Hsp_hit-to>([^\<]+)/){
      $blast{$id}{$hit}->{tend} = $1 if ($1 > $blast{$id}{$hit}->{tend});
    }
    elsif ($_ =~ m/<Hsp_query-from>([^\<]+)/){
      $blast{$id}{$hit}->{start} =$1;
    }
    elsif ($_ =~ m/<Hsp_query-to>([^<]+)/){
      $blast{$id}{$hit}->{end} = $1;
    }
    elsif ($_ =~ m/<Hsp_bit-score>([^\<]+)/){
      $blast{$id}{$hit}->{score} = $1;
    }
  }
  close IN;

  foreach my $contig (sort keys %blast){
 #   print "$contig\n";
    foreach my $h (sort {$a<=>$b} keys %{$blast{$contig}}){
      #print "$h\n";
      my $strand;
      if ($blast{$contig}{$h}->{tstart} > $blast{$contig}{$h}->{tend}){
        $strand = "-";
      }
      else {
        $strand = "+";
      }
      my @abf = split //, $ab[$i];
      my $fam = $abf[0] . $abf[1] . $abf[2];
      print "$contig\tresfinder\tblock\t$blast{$contig}{$h}->{start}\t$blast{$contig}{$h}->{end}\t$blast{$contig}{$h}->{score}\t$strand\t.\tID=$contig.$fam" . "$h\n";
      print "$contig\tresfinder\tblock1\t$blast{$contig}{$h}->{start}\t$blast{$contig}{$h}->{end}\t$blast{$contig}{$h}->{score}\t$strand\t.\tID=$contig.$fam" . "$h.match;Parent=$contig.$fam" . "$h;Note=$ab[$i]-resistance;Target=$blast{$contig}{$h}->{definition}\n";
    }
  }
  $i++;
}
