#!/usr/bin/perl

#Author: Jessica Gomez-Garrido, CNAG-CRG.
#Contact email: jessica.gomez@cnag.crg.eu
use strict;
use warnings;
use Getopt::Long;

my ($prokka, $rgi, $pref, $pep);

GetOptions(
           'prokka:s'        => \$prokka,
           'pep:s'           => \$pep,
           'rgi:s'           => \$rgi,
           'pref:s'          => \$pref
           );

my %proteins;
my @file = split /\n/, `FastaToTbl $pep`; 
my $len = scalar @file;
my $n = 0;
while ($n < $len){
  my @line = split /\s/, $file[$n];
  $proteins{$line[0]} = $line[1];
  $n++;
}


open RGI,"< $rgi" || die "cannot open input file $rgi";
my %rgi_annot; 
my %rgi_loc;
my $sample;
my %ori;
while (<RGI>){
  chomp;
  my @tmp = split //, $_;
  pop @tmp;
  $_ = join '', @tmp;
  my @line = split /\t/,$_;
  if ($line[0] eq "ORF_ID"){
    my $l = scalar @line;
    $rgi_loc{"header"}= "complete\tstart_type\trbs_motif\trbs_spacer\tgc_cont\tcut_off\tpass_bitscore\tbest_hit_bitscore\tbest_hit_aro\tbest_identities\taro\tmodel_type\tsnps_in_best_hit_aro\tother_snps\tdrug_class\tresistance_mechanism\tamr_gene_family\tpredicted_dna\tpredicted_protein\tcard_protein_sequence\tpercentage_length_of_reference_sequence";
  }
  else {
    my $id = $line[1];
    my @scaff = split /_/, $id;
    pop @scaff;
    my $contig = join "_", @scaff; 
    $sample = shift @scaff;
    my @info = split /\s/, $line[0];
    my $info = pop @info;
    $info .= ";";
    my $summary;
    if ($info =~ m/partial=([^;]+);start_type=([^;]+);rbs_motif=([^;]+);rbs_spacer=([^;]+);gc_cont=([^;]+)/){
      if ($1 eq "00"){
        $summary = "True\t";
      }
      else {
        $summary = "False\t";
      }
      $summary .= "$2\t$3\t$4\t$5";
    }
    $rgi_loc{$contig}{$line[2]}{$line[3]} = "$sample\t$contig\t$summary\t";
    $ori{$contig}{$line[2]}{$line[3]} = $line[4];
    my $l = 21;
    my $n = 5;
    while ($n < $l){
      $rgi_loc{$contig}{$line[2]}{$line[3]} .= "$line[$n]\t";
      $n++;
    }
    $rgi_annot{$contig}{$line[2]}{$line[3]} = $line[8];
    my $rest = shift @line;
  }
}
close RGI;

open PROKKA,"< $prokka" || die "cannot open input file $prokka";
open PROKKA_TBL, "> $pref.annot.tbl";
open JBROWSE, "> $pref.annot.4jb.gff3";
open RGI_TBL, "> $pref.rgi.tbl";
print RGI_TBL "gene\t" . $rgi_loc{"header"} . "\n";

print PROKKA_TBL "scaffold\tstart\tend\torientation\tgene\tgene_name\trgi\tec_number\tproduct\tinference\tjbrowse_link\tprotein_sequence\n";
my $total = keys %proteins; 
my @char = split //, $total;
my $i = scalar @char;
my %seen;
my $genes = 0;
while (<PROKKA>){
  chomp;
  next if ($_ =~ m/^#/);
  my @line = split /\t/, $_;
  next if (scalar @line != 9); 
  if ($line[2] eq "CDS"){
    $genes++;
    my $gene_id = sprintf("%0$i" . "d", $genes);
    my $newid = $pref . "g" . $gene_id;
    print PROKKA_TBL "$line[0]\t$line[3]\t$line[4]\t$line[6]\t$newid\t"; 
    my @jbline = @line;
    $jbline[1] = "CNAG";
    $jbline[2] = "gene";
    $jbline[8] = "ID=$newid";
    my $p = join "\t", @jbline;  
    print JBROWSE "$p\n";
    $jbline[2] = "mRNA";
    $jbline[8] = $line[8];
    $jbline[8] =~ s/ID=([^;]+);Parent=([^;]+)/ID=$newid.mrna;Parent=$newid/;    
    $line[8] .= ";";

    if ($line[8] =~ m/gene=([^;]+)/){
      print PROKKA_TBL "$1\t";
    }
    else {
      print PROKKA_TBL "\t";
    }
    if (exists $rgi_loc{$line[0]}{$line[3]}{$line[4]}){
      print PROKKA_TBL "Yes\t";
      $jbline[8] .= ";Description=$rgi_annot{$line[0]}{$line[3]}{$line[4]}";
      my @tbl= split /\t/, $rgi_loc{$line[0]}{$line[3]}{$line[4]};
      shift @tbl; 
      shift @tbl;
      my $tbl = join "\t", @tbl;
      print RGI_TBL "$newid\t$tbl\n";
      $seen{$line[0]}{$line[3]}{$line[4]}++;
      my @headers = split /\t/, $rgi_loc{"header"}; 
      my @info = split /\t/, $rgi_loc{$line[0]}{$line[3]}{$line[4]};
      my $h = 6;
      while ($h < scalar @info - 6){
        $jbline[8] .= ";$headers[$h]=$info[$h]";
        $h++;
      }
    }
    else {
      print PROKKA_TBL "\t";
    }
    if ($line[8] =~ m/eC_number=([^;]+)/){
      print PROKKA_TBL "$1\t";
    }
    else {
      print PROKKA_TBL "\t";
    }
    if ($line[8] =~ m/product=([^;]+)/){
      print PROKKA_TBL "$1\t";
    }
    else {
      print PROKKA_TBL "\t";
    }
    if ($line[8] =~ m/inference=([^;]+)/){
      print PROKKA_TBL "$1\t";
    }
    else {
      print PROKKA_TBL "\t";
    }
    print PROKKA_TBL "http://denovo.cnag.cat/genomes/cre/browse/?loc=$line[0]%3A$line[3]..$line[4]\t";
    if ($line[8] =~ m/ID=([^;]+)/){
      print PROKKA_TBL "$proteins{$1}\n";
    }
    $jbline[8] =~ s/locus_tag=([^;]+)/locus_tag=$newid/;
    $p = join "\t", @jbline;  
    print JBROWSE "$p\n";
    if (exists $rgi_loc{$line[0]}{$line[3]}{$line[4]}){
      $jbline[2] = "rgiCDS";
    }
    else {
      $jbline[2] = "CDS";
    }
    $jbline[8] =~ s/ID=([^;]+);Parent=([^;]+)/ID=$newid.cds;Parent=$newid.mrna/;

    $p = join "\t", @jbline;
    print JBROWSE "$p\n";
  }
}
close PROKKA;
foreach my $c (keys %rgi_loc){
  next if ($c eq "header");
  foreach my $s (keys %{$rgi_loc{$c}}){
    foreach my $e (keys %{$rgi_loc{$c}{$s}}){
      next if (exists $seen{$c}{$s}{$e});
      $genes++;
      my $gene_id = sprintf("%0$i" . "d", $genes);
      my $newid = $pref . "g" . $gene_id;
      my @line = split /\t/, $rgi_loc{$c}{$s}{$e};
      my $p = "$c\tCNAG\tgene\t$s\t$e\t.\t$ori{$c}{$s}{$e}\t.\tID=$newid";    
      print JBROWSE "$p\n";
      my @jbline = split /\t/, $p; 
      $jbline[2] = "mRNA";
      $jbline[8] .= ".mrna;Parent=$newid;locus_tag=$newid;Description=$rgi_annot{$c}{$s}{$e}"; 
      my @headers = split /\t/, $rgi_loc{"header"}; 
      my @info = split /\t/, $rgi_loc{$c}{$s}{$e};
      my $h = 6;
      while ($h < scalar @info - 6){
        $jbline[8] .= ";$headers[$h]=$info[$h]";
        $h++;
      }
      $p = join "\t",@jbline;
      print JBROWSE "$p\n";
      $jbline[2] = "rgiCDS";
      $jbline[8] =~ s/ID=([^;]+);Parent=([^;]+)/ID=$newid.cds;Parent=$newid.mrna/;
      $p = join "\t", @jbline;
      print JBROWSE "$p\n";
      my @sample = split /_/, $c;
      my @tbl= split /\t/, $rgi_loc{$c}{$s}{$e};
      shift @tbl; 
      shift @tbl;
      my $tbl = join "\t", @tbl;
      print RGI_TBL "$newid\t$tbl\n";
      print PROKKA_TBL "$c\t$s\t$e\t$ori{$c}{$s}{$e}\t$newid\t\tYes\t\t\t\thttp://denovo.cnag.cat/genomes/cre/browse/?loc=$c%3A$s..$e\t$line[20]\n"; 
    }
  }
}
close PROKKA_TBL;
close RGI_TBL;
close JBROWSE;


