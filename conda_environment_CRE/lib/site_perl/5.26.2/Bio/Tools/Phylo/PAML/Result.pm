package Bio::Tools::Phylo::PAML::Result;
$Bio::Tools::Phylo::PAML::Result::VERSION = '1.7.3';
use utf8;
use strict;
use warnings;

use base qw(Bio::Root::Root Bio::AnalysisResultI);

# ABSTRACT: A PAML result set object
# AUTHOR: Jason Stajich <jason@bioperl.org>
# AUTHOR: Aaron Mackey <amackey@virginia.edu>
# OWNER: Jason Stajich <jason@bioperl.org>
# OWNER: Aaron Mackey <amackey@virginia.edu>
# LICENSE: Perl_5

# AUTHOR: Albert Vilella <avilella@gmail.com>



sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($trees,$mlmat,$seqs,$ngmatrix,
      $codonpos,$codonfreq,$version,
      $model,$patterns, $stats,
      $aafreq, $aadistmat,
      $aamldistmat,
      $ntfreqs, $seqfile, $kappa_mat, $alpha_mat,
      $NSSitesresults,$input_params,$rst,$rst_persite,$rst_trees ) =
      $self->_rearrange([qw
                 (TREES MLMATRIX
                  SEQS NGMATRIX
                  CODONPOS CODONFREQ
                  VERSION MODEL PATTERNS
                  STATS AAFREQ AADISTMAT
                  AAMLDISTMAT
                  NTFREQ SEQFILE
                  KAPPA_DISTMAT
                  ALPHA_DISTMAT
                  NSSITESRESULTS
                  INPUT_PARAMS
                  RST RST_PERSITE RST_TREES)],
                @args);
  $self->reset_seqs;
  if( $trees ) {
      if(ref($trees) !~ /ARRAY/i ) {
      $self->warn("Must provide a valid array reference to initialize trees");
      } else {
      foreach my $t ( @$trees ) {
          $self->add_tree($t);
      }
      }
  }
  $self->{'_treeiterator'} = 0;

  if( $mlmat ) {
      if( ref($mlmat) !~ /ARRAY/i ) {
      $self->warn("Must provide a valid array reference to initialize MLmatrix");
      } else {
      $self->set_MLmatrix($mlmat);
      }
  }
  if( $seqs ) {
      if( ref($seqs) !~ /ARRAY/i ) {
      $self->warn("Must provide a valid array reference to initialize seqs");
      } else {
      foreach my $s ( @$seqs ) {
          $self->add_seq($s);
      }
      }
  }
  if( $ngmatrix ) {
      if( ref($ngmatrix) !~ /ARRAY/i ) {
      $self->warn("Must provide a valid array reference to initialize NGmatrix");
      } else {
      $self->set_NGmatrix($ngmatrix);
      }
  }
  if( $codonfreq ) {
      if( ref($codonfreq) =~ /ARRAY/i ) {
      $self->set_CodonFreqs($codonfreq);
      } else {
      $self->warn("Must provide a valid array reference to initialize codonfreq");
      }
  }

  if( $codonpos ) {
      if( ref($codonpos) !~ /ARRAY/i ) {
      $self->warn("Must provide a valid array reference to initialize codonpos");
      } else {
      $self->set_codon_pos_basefreq(@$codonpos);
      }
  }

  $self->version($version)   if defined $version;
  $self->seqfile($seqfile)   if defined $seqfile;
  $self->model($model)       if defined $model;
  if( defined $patterns ) {
      if( ref($patterns) =~ /HASH/i ) {
      $self->patterns($patterns);
      } else {
      $self->warn("Must provide a valid array reference to initialize patterns");
      }
  }

  $self->{'_aafreqs'} = {};
  if( $aafreq ) {
      if( ref($aafreq) =~ /HASH/i ) {
      $self->set_AAFreqs($aafreq);
      } else {
      $self->warn("Must provide a valid hash reference to initialize aafreq");
      }
  }
  if( $stats ) {
      if( ref($stats) =~ /HASH/i ) {
      while( my ($stat,$val) = each %$stats) {
          $self->add_stat($stat,$val);
      }
      } else {
      $self->warn("Must provide a valid hash reference initialize stats");
      }
  }
  $self->set_AADistMatrix($aadistmat) if defined $aadistmat;
  $self->set_AAMLDistMatrix($aamldistmat) if defined $aamldistmat;

  if( defined $NSSitesresults ) {
      if( ref($NSSitesresults) !~ /ARRAY/i ) {
      $self->warn("expected an arrayref for -NSSitesresults");
      } else {
      foreach my $m ( @$NSSitesresults ) {
          $self->add_NSSite_result($m);
      }
      }
  }

  $self->{'_ntfreqs'} = {};
  if( $ntfreqs ) {
      if( ref($ntfreqs) =~ /HASH/i ) {
      $self->set_NTFreqs($ntfreqs);
      } else {
      $self->warn("Must provide a valid hash reference to initialize ntfreq");
      }
  }

  if( $kappa_mat ) {
      $self->set_KappaMatrix($kappa_mat);
  }
  if( $alpha_mat ) {
      $self->set_AlphaMatrix($alpha_mat);
  }

  if( $input_params ) {
      if(  ref($input_params) !~ /HASH/i ) {
      $self->warn("Must provide a valid hash object for input_params\n");
      } else {
      while( my ($p,$v) = each %$input_params ) {
          $self->set_input_parameter($p,$v);
      }
      }

  }
  $self->reset_rst_seqs;
  if( $rst ) {
      if( ref($rst) =~ /ARRAY/i ) {
      for ( @$rst ) {
          $self->add_rst_seq($_);
      }
      } else {
      $self->warn("Need a valid array ref for -rst option\n");
      }
  }
  if( defined $rst_persite ) {
      $self->set_rst_persite($rst_persite);
  }
  $self->reset_rst_trees;
  if( $rst_trees ) {
      if( ref($rst_trees) =~ /ARRAY/i ) {
      for ( @$rst_trees ) {
          $self->add_rst_tree($_);
      }
      } else {
      $self->warn("Need a valid array ref for -rst_trees option\n");
      }
  }

  return $self;
}


sub next_tree{
   my ($self,@args) = @_;
   return $self->{'_trees'}->[$self->{'_treeiterator'}++] || undef;
}


sub get_trees{
   my ($self) = @_;
   return @{$self->{'_trees'} || []};
}


sub rewind_tree_iterator {
    shift->{'_treeiterator'} = 0;
}


sub add_tree{
   my ($self,$tree) = @_;
   if( $tree && ref($tree) && $tree->isa('Bio::Tree::TreeI') ) {
       push @{$self->{'_trees'}},$tree;
   }
   return scalar @{$self->{'_trees'}};
}



sub set_MLmatrix{
   my ($self,$mat) = @_;
   return unless ( defined $mat );
   if( ref($mat) !~ /ARRAY/i ) {
       $self->warn("Did not provide a valid 2D Array reference for set_MLmatrix");
       return;
   }
   $self->{'_mlmatrix'} = $mat;
}


sub get_MLmatrix{
   my ($self,@args) = @_;
   return $self->{'_mlmatrix'};
}


sub set_NGmatrix{
   my ($self,$mat) = @_;
   return unless ( defined $mat );
   if( ref($mat) !~ /ARRAY/i ) {
       $self->warn("Did not provide a valid 2D Array reference for set_NGmatrix");
       return;
   }
   $self->{'_ngmatrix'} = $mat;
}


sub get_NGmatrix{
   my ($self,@args) = @_;
   return $self->{'_ngmatrix'};
}



sub add_seq{
   my ($self,$seq) = @_;
   if( $seq ) {
       unless( $seq->isa("Bio::PrimarySeqI") ) {
       $self->warn("Must provide a valid Bio::PrimarySeqI to add_seq");
       return;
       }
       push @{$self->{'_seqs'}},$seq;
   }

}


sub reset_seqs{
   my ($self) = @_;
   $self->{'_seqs'} = [];
}


sub get_seqs{
   my ($self) = @_;
   return @{$self->{'_seqs'}};
}


sub set_codon_pos_basefreq {
    my ($self,@codonpos) = @_;
    if( scalar @codonpos != 3 ) {
    $self->warn("invalid array to set_codon_pos_basefreq, must be an array of length 3");
    return;
    }
    foreach my $pos ( @codonpos ) {
    if( ref($pos) !~ /HASH/i ||
        ! exists $pos->{'A'} ) {
        $self->warn("invalid array to set_codon_pos_basefreq, must be an array with hashreferences keyed on DNA bases, C,A,G,T");
    }
    }
    $self->{'_codonposbasefreq'} = [@codonpos];
}


sub get_codon_pos_basefreq{
   my ($self) = @_;
   return @{$self->{'_codonposbasefreq'}};
}


sub version{
   my $self = shift;
   $self->{'_version'} = shift if @_;
   return $self->{'_version'};
}


sub seqfile{
   my $self = shift;
   $self->{'_seqfile'} = shift if @_;
   return $self->{'_seqfile'};
}


sub model{
    my $self = shift;

    return $self->{'_model'} = shift if @_;
    return $self->{'_model'};
}



sub patterns{
    my $self = shift;
    return $self->{'_patterns'} = shift if @_;
    return $self->{'_patterns'};
}


sub set_AAFreqs{
   my ($self,$aafreqs) = @_;

   if( $aafreqs && ref($aafreqs) =~ /HASH/i ) {
       foreach my $seqname ( keys %{$aafreqs} ) {
       $self->{'_aafreqs'}->{$seqname} = $aafreqs->{$seqname};
       }
   }
}


sub get_AAFreqs{
   my ($self,$seqname) = @_;
   if( $seqname ) {
       return $self->{'_aafreqs'}->{$seqname} || {};
   } else {
       return $self->{'_aafreqs'};
   }
}


sub set_NTFreqs{
   my ($self,$freqs) = @_;

   if( $freqs && ref($freqs) =~ /HASH/i ) {
       foreach my $seqname ( keys %{$freqs} ) {
       $self->{'_ntfreqs'}->{$seqname} = $freqs->{$seqname};
       }
   }
}


sub get_NTFreqs{
   my ($self,$seqname) = @_;
   if( $seqname ) {
       return $self->{'_ntfreqs'}->{$seqname} || {};
   } else {
       return $self->{'_ntfreqs'};
   }
}


sub add_stat{
   my ($self,$stat,$value) = @_;
   return if( ! defined $stat || !defined $value );
   $self->{'_stats'}->{$stat} = $value;
   return;
}


sub get_stat{
   my ($self,$statname) = @_;
   return $self->{'_stats'}->{$statname};
}


sub get_stat_names{
   my ($self) = @_;
   return keys %{$self->{'_stats'} || {}};
}


sub get_AADistMatrix{
    my $self = shift;
    return $self->{'_AADistMatix'};
}


sub set_AADistMatrix{
   my ($self,$d) = @_;
   if( ! $d ||
       ! ref($d) ||
       ! $d->isa('Bio::Matrix::PhylipDist') ) {
       $self->warn("Must provide a valid Bio::Matrix::MatrixI for set_AADistMatrix");
   }
   $self->{'_AADistMatix'} = $d;
   return;
}


sub get_AAMLDistMatrix{
    my $self = shift;
    return $self->{'_AAMLDistMatix'};
}


sub set_AAMLDistMatrix{
   my ($self,$d) = @_;
   if( ! $d ||
       ! ref($d) ||
       ! $d->isa('Bio::Matrix::PhylipDist') ) {
       $self->warn("Must provide a valid Bio::Matrix::MatrixI for set_AAMLDistMatrix");
   }
   $self->{'_AAMLDistMatix'} = $d;
   return;
}


sub add_NSSite_result{
   my ($self,$model) = @_;
   if( defined $model ) {
       push @{$self->{'_nssiteresult'}}, $model;
   }
   return scalar @{$self->{'_nssiteresult'}};
}


sub get_NSSite_results{
   my ($self) = @_;
   return @{$self->{'_nssiteresult'} || []};
}


sub set_CodonFreqs{
    my $self = shift;

    return $self->{'_codonfreqs'} = shift if @_;
    return $self->{'_codonfreqs'};
}


sub get_CodonFreqs{
   my ($self) = @_;
   return @{$self->{'_codonfreqs'} || []};
}



sub get_KappaMatrix{
    my $self = shift;
    return $self->{'_KappaMatix'};
}


sub set_KappaMatrix{
   my ($self,$d) = @_;
   if( ! $d ||
       ! ref($d) ||
       ! $d->isa('Bio::Matrix::PhylipDist') ) {
       $self->warn("Must provide a valid Bio::Matrix::MatrixI for set_NTDistMatrix");
   }
   $self->{'_KappaMatix'} = $d;
   return;
}



sub get_AlphaMatrix{
    my $self = shift;
    return $self->{'_AlphaMatix'};
}


sub set_AlphaMatrix{
   my ($self,$d) = @_;
   if( ! $d ||
       ! ref($d) ||
       ! $d->isa('Bio::Matrix::PhylipDist') ) {
       $self->warn("Must provide a valid Bio::Matrix::MatrixI for set_NTDistMatrix");
   }
   $self->{'_AlphaMatix'} = $d;
   return;
}


sub set_input_parameter{
   my ($self,$p,$v) = @_;
   return unless defined $p;
   $self->{'_input_parameters'}->{$p} = $v;
}


sub get_input_parameters{
   my ($self) = @_;
   return %{$self->{'_input_parameters'} || {}};
}


sub reset_input_parameters{
   my ($self) = @_;
   $self->{'_input_parameters'} = {};
}


sub add_rst_seq{
   my ($self,$seq) = @_;
   if( $seq ) {
       unless( $seq->isa("Bio::PrimarySeqI") ) {
       $self->warn("Must provide a valid Bio::PrimarySeqI to add_rst_seq");
       return;
       }
       push @{$self->{'_rstseqs'}},$seq;
   }

}


sub reset_rst_seqs{
   my ($self) = @_;
   $self->{'_rstseqs'} = [];
}


sub get_rst_seqs{
   my ($self) = @_;
   return @{$self->{'_rstseqs'} || []};
}



sub add_rst_tree{
   my ($self,$tree) = @_;
   if( $tree ) {
       unless( $tree->isa("Bio::Tree::TreeI") ) {
       $self->warn("Must provide a valid Bio::Tree::TreeI to add_rst_tree not $tree");
       return;
       }
       push @{$self->{'_rsttrees'}},$tree;
   }
}


sub reset_rst_trees{
   my ($self) = @_;
   $self->{'_rsttrees'} = [];
}


sub get_rst_trees{
   my ($self) = @_;
   return @{$self->{'_rsttrees'} || []};
}


sub set_rst_persite{
    my $self = shift;

    return $self->{'_rstpersite'} = shift if @_;
    return $self->{'_rstpersite'};
}


sub get_rst_persite{
   my ($self) = @_;
   return $self->{'_rstpersite'} || [];
}



1;

__END__

=pod

=encoding UTF-8

=head1 NAME

Bio::Tools::Phylo::PAML::Result - A PAML result set object

=head1 VERSION

version 1.7.3

=head1 SYNOPSIS

  # see Bio::Tools::Phylo::PAML for example usage
  use Bio::Tools::Phylo::PAML;
  my $parser = Bio::Tools::Phylo::PAML->new
    (-file => "./results/mlc", -dir => "./results/");

  # get the first/next result; a Bio::Tools::Phylo::PAML::Result object,
  # which isa Bio::SeqAnalysisResultI object.
  my $result = $parser->next_result();

  my @seqs         = $result->get_seqs;
  my %input_params = $result->get_input_parameters;
  my @basfreq      = $result->get_codon_pos_basefreq;
  my $MLmatrix     = $result->get_MLmatrix; # get MaxLikelihood Matrix
  my $NGmatrix     = $result->get_NGmatrix; # get Nei-Gojoburi Matrix


  # for AAML runs
  my $AAmatrix   = $result->get_AADistMatrix;
  my $AAMLmatrix   = $result->get_AAMLDistMatrix;

  # if -dir contains an rst file get list of
  # Bio::PrimarySeq ancestral state reconstructions of the sequences
  my @rsts          = $result->get_rst_seqs;


  # if you want to print the changes on the tree
  # this will print out the
  # anc_aa       => ANCESTRAL AMINO ACID
  # anc_prob     => ANCESTRAL AA PROBABILITY
  # derived_aa   => DERIVED AA
  # derived_prob => DERIVE AA PROBABILITY (where appropriate - NA for extant/tip taxas)
  # site         => which codon site this in the alignment
    @trees = $result->get_rst_trees;
    for my $t ( @trees ) {
    for my $node ( $t->get_nodes ) {
        next unless $node->ancestor; # skip root node
        my @changes = $node->get_tag_values('changes');
        my $chgstr = '';
        for my $c ( @changes ) {
        for my $k ( sort keys %$c ) {
            $chgstr .= "$k => $c->{$k} ";
        }
        $chgstr .= "\n\t";
        }

        printf "node:%s n=%s s=%s\n\t%s\n",
        $node->id,
        $node->get_tag_values('n'),
        $node->get_tag_values('s'),
        $chgstr;
    }
    }

  # Persite probabilities
  my $persite = $result->get_rst_persite;
  # let's score site 1
  $site = $persite->[2];
  # so site 2, node 2 (extant node, node 2)
  print $site->[2]->{'codon'}, ' ',$site->[2]->{'aa'},"\n";
  # site 2, node 3
  print $site->[3]->{'codon'}, ' ',$site->[3]->{'aa'}, "\n";

  # ancestral node 9, codon, aa, marginal probabilities; Yang95 is listed as
  #  (eqn. 4 in Yang et al. 1995 Genetics 141:1641-1650) in PAML rst file.
  print $site->[9]->{'codon'}, ' ',$site->[9]->{'aa'}, ' ', $site->[9]->{'prob'}, ' ',
        $site->[9]->{'Yang95_aa'},' ', $site->[9]->{'Yang95_aa_prob'},"\n";

=head1 DESCRIPTION

This is a container object for PAML Results.

=head1 METHODS

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Phylo::PAML::Result->new(%data);
 Function: Builds a new Bio::Tools::Phylo::PAML::Result object
 Returns : Bio::Tools::Phylo::PAML::Result
 Args    : -trees     => array reference of Bio::Tree::TreeI objects
           -MLmatrix  => ML matrix
           -seqs      => array reference of Bio::PrimarySeqI objects
           -codonpos  => array reference of codon positions
           -codonfreq => array reference of codon frequencies
           -version   => version string
           -model     => model string
           -patterns  => hashref with the fields '-patterns', '-ns', '-ls'
           -stats     => array ref of misc stats   (optional)
           -aafreq    => Hashref of AA frequencies (only for AAML)
           -aadistmat => Bio::Matrix::PhylipDist   (only for AAML)
           -aamldistmat => Bio::Matrix::PhylipDist   (only for pairwise AAML)
           -ntfreq    => array ref of NT frequencies (only for BASEML)
           -seqfile    => seqfile used
           -kappa_mat => Bio::Matrix::PhylipDist of kappa values (only for BASEML)
           -alpha_mat => Bio::Matrix::PhylipDist of alpha values (only for BASEML)
           -NSSitesresult => arrayref of PAML::ModelResult
           -input_params  => input params from .ctl file
           -rst       => array reference of Bio::PrimarySeqI objects
                         of ancestral state reconstruction
           -rst_persite=> arrayref of persite data, this is a complicated set of AoH
           -rst_trees  => rst trees with changes coded on the tree

See Also: L<Bio::Tree::TreeI>, L<Bio::PrimarySeqI>, L<Bio::Matrix::PhylipDist>, L<Bio::Tools::Phylo::PAML>

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $factory->next_tree;
 Function: Get the next tree from the factory
 Returns : L<Bio::Tree::TreeI>
 Args    : none

=head2 get_trees

 Title   : get_trees
 Usage   : my @trees = $result->get_trees;
 Function: Get all the parsed trees as an array
 Returns : Array of trees
 Args    : none

=head2 rewind_tree_iterator

 Title   : rewind_tree_iterator
 Usage   : $result->rewind_tree_iterator()
 Function: Rewinds the tree iterator so that next_tree can be
           called again from the beginning
 Returns : none
 Args    : none

=head2 add_tree

 Title   : add_tree
 Usage   : $result->add_tree($tree);
 Function: Adds a tree
 Returns : integer which is the number of trees stored
 Args    : L<Bio::Tree::TreeI>

=head2 set_MLmatrix

 Title   : set_MLmatrix
 Usage   : $result->set_MLmatrix($mat)
 Function: Set the ML Matrix
 Returns : none
 Args    : Arrayref to MLmatrix (must be arrayref to 2D matrix whic is
       lower triangle pairwise)

=head2 get_MLmatrix

 Title   : get_MLmatrix
 Usage   : my $mat = $result->get_MLmatrix()
 Function: Get the ML matrix
 Returns : 2D Array reference
 Args    : none

=head2 set_NGmatrix

 Title   : set_NGmatrix
 Usage   : $result->set_NGmatrix($mat)
 Function: Set the Nei & Gojobori Matrix
 Returns : none
 Args    : Arrayref to NGmatrix (must be arrayref to 2D matrix whic is
       lower triangle pairwise)

=head2 get_NGmatrix

 Title   : get_NGmatrix
 Usage   : my $mat = $result->get_NGmatrix()
 Function: Get the Nei & Gojobori matrix
 Returns : 2D Array reference
 Args    : none

=head2 add_seq

 Title   : add_seq
 Usage   : $obj->add_seq($seq)
 Function: Add a Bio::PrimarySeq to the Result
 Returns : none
 Args    : Bio::PrimarySeqI
See also : L<Bio::PrimarySeqI>

=head2 reset_seqs

 Title   : reset_seqs
 Usage   : $result->reset_seqs
 Function: Reset the OTU seqs stored
 Returns : none
 Args    : none

=head2 get_seqs

 Title   : get_seqs
 Usage   : my @otus = $result->get_seqs
 Function: Get the seqs Bio::PrimarySeq (OTU = Operational Taxonomic Unit)
 Returns : Array of Bio::PrimarySeq
 Args    : None
See also : L<Bio::PrimarySeq>

=head2 set_codon_pos_basefreq

 Title   : set_codon_pos_basefreq
 Usage   : $result->set_codon_pos_basefreq(@freqs)
 Function: Set the codon position base frequencies
 Returns : none
 Args    : Array of length 3 where each slot has a hashref
           keyed on DNA base

=head2 get_codon_pos_basefreq

 Title   : get_codon_pos_basefreq
 Usage   : my @basepos = $result->get_codon_pos_basefreq;
 Function: Get the codon position base frequencies
 Returns : Array of length 3 (each codon position), each
           slot is a hashref keyed on DNA bases, the values are
           the frequency of the base at that position for all sequences
 Args    : none
 Note    : The array starts at 0 so position '1' is in position '0'
           of the array

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: Get/Set version
 Returns : value of version
 Args    : newvalue (optional)

=head2 seqfile

 Title   : seqfile
 Usage   : $obj->seqfile($newval)
 Function: Get/Set seqfile
 Returns : value of seqfile
 Args    : newvalue (optional)

=head2 model

 Title   : model
 Usage   : $obj->model($newval)
 Function: Get/Set model
 Returns : value of model
 Args    : on set, new value (a scalar or undef, optional)

=head2 patterns

 Title   : patterns
 Usage   : $obj->patterns($newval)
 Function: Get/Set Patterns hash
 Returns : Hashref of pattern data
 Args    : [optional] Hashref of patterns
         : The hashref is typically
         : { -patterns => \@arrayref
         :   -ns       => $ns
         :   -ls       => $ls
         : }

=head2 set_AAFreqs

 Title   : set_AAFreqs
 Usage   : $result->set_AAFreqs(\%aafreqs);
 Function: Get/Set AA freqs
 Returns : none
 Args    : Hashref, keys are the sequence names, each points to a hashref
           which in turn has keys which are the amino acids

=head2 get_AAFreqs

 Title   : get_AAFreqs
 Usage   : my %all_aa_freqs = $result->get_AAFreqs()
            OR
           my %seq_aa_freqs = $result->get_AAFreqs($seqname)
 Function: Get the AA freqs, either for every sequence or just
           for a specific sequence
           The average aa freqs for the entire set are also available
           for the sequence named 'Average'
 Returns : Hashref
 Args    : (optional) sequence name to retrieve aa freqs for

=head2 set_NTFreqs

 Title   : set_NTFreqs
 Usage   : $result->set_NTFreqs(\%aafreqs);
 Function: Get/Set NT freqs
 Returns : none
 Args    : Hashref, keys are the sequence names, each points to a hashref
           which in turn has keys which are the amino acids

=head2 get_NTFreqs

 Title   : get_NTFreqs
 Usage   : my %all_nt_freqs = $result->get_NTFreqs()
            OR
           my %seq_nt_freqs = $result->get_NTFreqs($seqname)
 Function: Get the NT freqs, either for every sequence or just
           for a specific sequence
           The average nt freqs for the entire set are also available
           for the sequence named 'Average'
 Returns : Hashref
 Args    : (optional) sequence name to retrieve nt freqs for

=head2 add_stat

 Title   : add_stat
 Usage   : $result->add_stat($stat,$value);
 Function: Add some misc stat valuess (key/value pairs)
 Returns : none
 Args    : $stat  stat name
           $value stat value

=head2 get_stat

 Title   : get_stat
 Usage   : my $value = $result->get_stat($name);
 Function: Get the value for a stat of a given name
 Returns : scalar value
 Args    : name of the stat

=head2 get_stat_names

 Title   : get_stat_names
 Usage   : my @names = $result->get_stat_names;
 Function: Get the stat names stored for the result
 Returns : array of names
 Args    : none

=head2 get_AADistMatrix

 Title   : get_AADistMatrix
 Usage   : my $mat = $obj->get_AADistMatrix()
 Function: Get AADistance Matrix
 Returns : value of AADistMatrix (Bio::Matrix::PhylipDist)
 Args    : none

=head2 set_AADistMatrix

 Title   : set_AADistMatrix
 Usage   : $obj->set_AADistMatrix($mat);
 Function: Set the AADistrance Matrix (Bio::Matrix::PhylipDist)
 Returns : none
 Args    : AADistrance Matrix (Bio::Matrix::PhylipDist)

=head2 get_AAMLDistMatrix

 Title   : get_AAMLDistMatrix
 Usage   : my $mat = $obj->get_AAMLDistMatrix()
 Function: Get AAMLDistance Matrix
 Returns : value of AAMLDistMatrix (Bio::Matrix::PhylipDist)
 Args    : none

=head2 set_AAMLDistMatrix

 Title   : set_AAMLDistMatrix
 Usage   : $obj->set_AAMLDistMatrix($mat);
 Function: Set the AA ML Distrance Matrix (Bio::Matrix::PhylipDist)
 Returns : none
 Args    : AAMLDistrance Matrix (Bio::Matrix::PhylipDist)

=head2 add_NSSite_result

 Title   : add_NSSite_result
 Usage   : $result->add_NSSite_result($model)
 Function: Add a NSsite result (PAML::ModelResult)
 Returns : none
 Args    : Bio::Tools::Phylo::PAML::ModelResult

=head2 get_NSSite_results

 Title   : get_NSSite_results
 Usage   : my @results = @{$self->get_NSSite_results};
 Function: Get the reference to the array of NSSite_results
 Returns : Array of PAML::ModelResult results
 Args    : none

=head2 set_CodonFreqs

 Title   : set_CodonFreqs
 Usage   : $obj->set_CodonFreqs($newval)
 Function: Get/Set the Codon Frequence table
 Returns : value of set_CodonFreqs (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=head2 get_CodonFreqs

 Title   : get_CodonFreqs
 Usage   : my @codon_freqs = $result->get_CodonFreqs()
 Function: Get the Codon freqs
 Returns : Array
 Args    : none

=head2 get_KappaMatrix

 Title   : get_KappaMatrix
 Usage   : my $mat = $obj->get_KappaMatrix()
 Function: Get KappaDistance Matrix
 Returns : value of KappaMatrix (Bio::Matrix::PhylipDist)
 Args    : none

=head2 set_KappaMatrix

 Title   : set_KappaMatrix
 Usage   : $obj->set_KappaMatrix($mat);
 Function: Set the KappaDistrance Matrix (Bio::Matrix::PhylipDist)
 Returns : none
 Args    : KappaDistrance Matrix (Bio::Matrix::PhylipDist)

=head2 get_AlphaMatrix

 Title   : get_AlphaMatrix
 Usage   : my $mat = $obj->get_AlphaMatrix()
 Function: Get AlphaDistance Matrix
 Returns : value of AlphaMatrix (Bio::Matrix::PhylipDist)
 Args    : none

=head2 set_AlphaMatrix

 Title   : set_AlphaMatrix
 Usage   : $obj->set_AlphaMatrix($mat);
 Function: Set the AlphaDistrance Matrix (Bio::Matrix::PhylipDist)
 Returns : none
 Args    : AlphaDistrance Matrix (Bio::Matrix::PhylipDist)

=head2 set_input_parameter

 Title   : set_input_parameter
 Usage   : $obj->set_input_parameter($p,$vl);
 Function: Set an Input Parameter
 Returns : none
 Args    : $parameter and $value

=head2 get_input_parameters

 Title   : get_input_parameters
 Usage   : $obj->get_input_parameters;
 Function: Get Input Parameters
 Returns : Hash of key/value pairs
 Args    : none

=head2 reset_input_parameters

 Title   : reset_input_parameters
 Usage   : $obj->reset_input_parameters;
 Function: Reset the Input Parameters hash
 Returns : none
 Args    : none

=head2 add_rst_seq

 Title   : add_rst_seq
 Usage   : $obj->add_rst_seq($seq)
 Function: Add a Bio::PrimarySeq to the RST Result
 Returns : none
 Args    : Bio::PrimarySeqI
See also : L<Bio::PrimarySeqI>

=head2 reset_rst_seqs

 Title   : reset_rst_seqs
 Usage   : $result->reset_rst_seqs
 Function: Reset the RST seqs stored
 Returns : none
 Args    : none

=head2 get_rst_seqs

 Title   : get_rst_seqs
 Usage   : my @otus = $result->get_rst_seqs
 Function: Get the seqs Bio::PrimarySeq
 Returns : Array of Bio::PrimarySeqI objects
 Args    : None
See also : L<Bio::PrimarySeq>

=head2 add_rst_tree

 Title   : add_rst_tree
 Usage   : $obj->add_rst_tree($tree)
 Function: Add a Bio::Tree::TreeI to the RST Result
 Returns : none
 Args    : Bio::Tree::TreeI
See also : L<Bio::Tree::TreeI>

=head2 reset_rst_trees

 Title   : reset_rst_trees
 Usage   : $result->reset_rst_trees
 Function: Reset the RST trees stored
 Returns : none
 Args    : none

=head2 get_rst_trees

 Title   : get_rst_trees
 Usage   : my @otus = $result->get_rst_trees
 Function: Get the trees Bio::Tree::TreeI
 Returns : Array of Bio::Tree::TreeI objects
 Args    : None
See also : L<Bio::Tree::TreeI>

=head2 set_rst_persite

 Title   : set_rst_persite
 Usage   : $obj->set_rst_persite($newval)
 Function: Get/Set the per-site RST values
 Returns : value of set_rst_persite (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=head2 get_rst_persite

 Title   : get_rst_persite
 Usage   : my @rst_persite = @{$result->get_rst_persite()}
 Function: Get the per-site RST values
 Returns : Array
 Args    : none

=head1 Reconstructed Ancestral State relevant options

=head1 FEEDBACK

=head2 Mailing lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/Support.html    - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:
I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bio-tools-phylo-paml/issues

=head1 AUTHORS

Jason Stajich <jason@bioperl.org>

Aaron Mackey <amackey@virginia.edu>

Albert Vilella <avilella@gmail.com>

=head1 COPYRIGHT

This software is copyright (c) by Jason Stajich <jason@bioperl.org>, and by Aaron Mackey <amackey@virginia.edu>.

This software is available under the same terms as the perl 5 programming language system itself.

=cut
