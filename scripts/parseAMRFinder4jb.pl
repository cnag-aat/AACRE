#!/usr/bin/perl

#Author: Jessica Gomez-Garrido, CNAG-CRG.
#Contact email: jessica.gomez@cnag.crg.eu
use strict;
use warnings;
use Getopt::Long;

my ($annot, $amrf);

GetOptions(
           'annot:s'        => \$annot,
           'amrf:s'           => \$amrf
           );

open AMRF,"< $amrf" || die "cannot open input file $amrf";
my %amrf_annot;
while (<AMRF>){
  chomp;
  my @line = split /\t/,$_;
  if ($line[0] eq "gene"){
    my $l = scalar @line;
    $amrf_annot{"header"}= $_;
  }
  else {
    my $id = $line[0];
    $amrf_annot{$id} = $_;
  }
}
close AMRF; 

open ANNOT,"< $annot" || die "cannot open input file $annot";
my $gene;
while (<ANNOT>){
  chomp;
  my @line = split /\t/,$_;
  if ($line[2] eq 'gene' && $line[8] =~ m/ID=/){
    $gene = $';
    print "$_\n" if (exists $amrf_annot{$gene}); 
  }
  elsif ($line[2] =~ m/CDS/ && exists $amrf_annot{$gene}){
    $line[2]="amrfCDS";
    my @headers = split /\t/, $amrf_annot{"header"};
    my @info = split /\t/, $amrf_annot{$gene};
    $line[8] = "ID=$gene.cds;Parent=$gene.mrna;locus_tag=$gene;Name=$info[1]";
    my $h = 2;
    while ($h < scalar @info ){
      $line[8] .= ";$headers[$h]=$info[$h]";
      $h++;
    }
    my $p = join "\t", @line;
    print "$p\n";
  }
  elsif ($line[2] eq 'mRNA' && exists $amrf_annot{$gene}){
    my @headers = split /\t/, $amrf_annot{"header"};
    my @info = split /\t/, $amrf_annot{$gene};
    $line[8] = "ID=$gene.mrna;Parent=$gene;locus_tag=$gene;Name=$info[1]";
    my $h = 2;
    while ($h < scalar @info ){
      $line[8] .= ";$headers[$h]=$info[$h]";
      $h++;
    }
    my $p = join "\t", @line;
    print "$p\n";
  }
}
close ANNOT;
