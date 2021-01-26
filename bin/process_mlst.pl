#!/usr/bin/perl

#Author: Jessica Gomez-Garrido, CNAG-CRG.
#Contact email: jessica.gomez@cnag.crg.eu
use strict;
use warnings;



my $assembly = $ARGV[0];

while (<STDIN>){
  chomp;
  my @line = split /\t/, $_;
  print "assembly\tPubMLST\tST\talleles\n";
  if ($line[1] eq "-"){
    print "$assembly\t\t\t\n";
  }
  else {
    print "$assembly\t$line[1]\t$line[2]\t";
    my $i = 3;
    while ($i < scalar @line - 1){
      print "$line[$i],";
      $i++;
    }
    print "$line[$i]\n";
  }
}
