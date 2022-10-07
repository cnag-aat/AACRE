#! /usr/bin/perl

# Program to sort reversly and print out the bracken abundances from STDIN. Threshold of minimal abundance has been set to 5% (0.05) by default
# USAGE:   cat  out/AI3013.NB21.461.4196AC.FAK60492.CRE_05.abundance.bracken |  revsort_bracken.pl
#          default threshold for minimum abundance otherwise set it as: 
#          cat bracken abundance_file \| revsort_bracken.pl \-min_abun value\n"
# OUTPUT FORMAT is ONESTRING WITH all species:abundance;....


# Fernando Cruz, 03-07-2019

use strict;
use Getopt::Long;
use File::Basename;

# define variables
my %abundance=(); # hash list
my $value; # value
my $species; # key
my @line;# array for parsing elements

my $min=0.05;# default abundance threshold
my $file=""; 

# set option for abundance threshold 
GetOptions(
           'f:s'  => \$file,
           'min:s'   => \$min # it's necessary to define the kind of variable :s means string. but perl automatically turns this into a number if its required
);
#default threshold for minimum abundance otherwise set it as: bracken-top-hits.pl \-f \-min_abun value\n";

if ($file eq ""){
       die print "provide abundance filename! bracken-top-hits.pl \-f filename \-min value\n"; 
}

my ($basename,$path,$ext) = fileparse($file,qw(\.bracken \.abundance));

# reading FILE
open (IN, "$file") or die "print cannot open $file\n";

while(<IN>)
{

 #store the hash
 chomp;
 next if ($_ =~ /fraction/); # skip header
 @line=split /\t/, $_;
 $species=$line[0];    
 $value=$line[6];   
 next if ($value <= $min);
 $abundance{$species}=$value;

}# END OF FILE 

close (IN) or die "print cannot open $file\n";

#Print base filename
    print "$basename\t";    
#Sort the keys of the hash according to the values
foreach $species (reverse sort { $abundance{$a} <=> $abundance{$b} } keys %abundance) {
    print "$species\:$abundance{$species}\;";
}
#separate each sample entry
print "\n";

exit;






