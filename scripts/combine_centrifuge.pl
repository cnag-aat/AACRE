#!/usr/bin/perl

#Author: Jessica Gomez-Garrido, CNAG-CRG.
#Contact email: jessica.gomez@cnag.crg.eu
use strict;
use warnings;

my $cout = $ARGV[0];
my $report = $ARGV[1];
my $species = $ARGV[2];
my $scaffolds = $ARGV[3];

#print STDERR "$report\n";
my @sample = split /\./, $cout;

open SPECIES, "<$species";
my %lims_sp;
while (<SPECIES>){
  chomp;
  my @line = split /\t/, $_;
  $lims_sp{$sample[0]} = $line[2];
}
close SPECIES;

open REPORT, "<$report";
my %ids;
while (<REPORT>){
  chomp;
  next if ($_ =~m/taxID/);
  my @line = split /\t/, $_;
  my @name = split /\s/, $line[0];
  next if ($line[2] =~ m/genus|family|order/);
  $ids{$line[1]}="$name[0] $name[1]";
}
close REPORT;

open SCAFF, "<$scaffolds";
my %circular;
while (<SCAFF>){
  chomp;
  my @line = split /\t/, $_;
  $circular{$line[1]} = $line[4];  
}
close SCAFF;

open OUT, "<$cout";
while (<OUT>){
  chomp;
  next if ($_=~ m/readID/);
  my @line = split /\t/, $_;
  if (!exists $ids{$line[2]}){
    print "$line[0]\t$circular{$line[0]}\t-\t-\t-\t$lims_sp{$sample[0]}\tFalse\tFalse\n";
  }
  else {
    my @species = split " ", $lims_sp{$sample[0]};
    my @hit = split " ", $ids{$line[2]};
    my $len_rel = $line[5]/$line[6];
    print "$line[0]\t$circular{$line[0]}\t$ids{$line[2]}\t$line[1]\t$len_rel\t$lims_sp{$sample[0]}\t";
    if ($species[0] eq $hit[0]){
      print "True\t";
    }
    else {
      print "False\t";
    }
    if ($species[0] eq $hit[0] && $species[1] eq $hit[1]){
      print "True\n";
    }
    else {
      print "False\n";
    }
  }
}
close OUT; 
