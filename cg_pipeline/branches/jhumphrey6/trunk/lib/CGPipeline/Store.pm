# $Id: $

=head1 NAME

CGPipeline::Store - a library for CRUD operations on a Bio::SeqFeature::Store feature database

=head1 AUTHOR

Jay Humphrey <jhumphrey6@gmail.com>

=cut

=head1 SYNOPSIS

=head1 TODO

finish perldoc...

=cut

package CGPipeline::Store; 
use base qw/CGPipeline/;
use Bio::SeqIO;
use Data::Dumper;
use Bio::DB::SeqFeature::Store;
use Bio::Tools::GFF;
use Bio::DB::GFF::Util::Rearrange 'rearrange';

sub new{
  my $class = shift;
  my %params = @_;
  die "new or existing sqlite database file must be defined\n" unless defined $params{'file'};
  my $object = bless({
    file=>$params{'file'},
    verbose=>$params{'verbose'} || 0,
    volatile=>$params{'volatile'} || 0,
    database=>undef
  },$class);
  $object->_connect;
  $object->load_config;#load the default config file or an alternative
  return $object;
}
sub logmsg{
  my $self = shift;
#if(! $self->{verbose} ){return;}
  my $msg = shift;
  printf STDERR "$msg";
}
sub db{
  my $self = shift;
  return $self->{'database'};
}
sub _connect{
  my $self=shift;
  my $is_new = 0;
  my $sqlitedb = $self->{'file'};
  if( ! -e $sqlitedb){$is_new=1;logmsg("creating new $sqlitedb\n");}
  $self->{'database'} = Bio::DB::SeqFeature::Store->new( -adaptor=>'DBI::SQLite',-dsn=>"$sqlitedb", -create=>$is_new, -autoindex=>1, -write=>1);
}
sub sane_id{
  my ($self,$seq_id)=@_;
  my $ok = 1;
  if ($seq_id =~ /\s+/){
    $ok = 0;
  }
  return $ok;
}
 
sub add_seq{
  my ($self,$seq)=@_;
  if( $self->seq_exists($seq->{'seq_id'})){
    $self->logmsg(sprintf("Warning:%s already exists\n",$seq->{'seq_id'}));
    return 0 unless $self->{volatile};
    $self->logmsg("Volatile mode!\n");
  }
  if( ! $self->sane_id($seq->{'seq_id'})){
    $self->logmsg(sprintf("Failed:seq_id %s is unsafe\n",$seq->{'seq_id'}));
    return 0;
  }
  return $self->db->insert_sequence($seq->{'seq_id'},$seq->{'seq'},0);
}
sub seq_exists{ # in certain cases it could read a cached list of ids
  my ($self,$seqid)=@_;
  my @ids = $self->db->seq_ids();
  return scalar grep(/^$seqid$/,@ids);
}
sub add_feature{
  my $self=shift; #feature should be a Bio::SeqFeatureI
  my $feat=shift;
  my $seq_id = $feat->seq_id;
  my ($locus_id) = $feat->get_tag_values('locus_id');
  if(!$self->seq_exists($seq_id)){
    $self->logmsg("Warning:sequence $seq_id does not exist\n");
  }
  my ($feature) = $self->db->features(-attributes=>{'locus_id'=>$locus_id});#,-seq_id=>$seq_id);
  if($feature){
    if($self->{volatile}){
      $self->logmsg("Replacing feature $seq_id:$locus_id\n");
      $self->db->delete(($feature));
    }
    else{logmsg ("skipping $locus_id (exists)\n");return $feature;}
  }
  return $self->db->store($feat);
}
sub write_gff{
  my ($self,$gfffile)=@_;
  my $gff = Bio::Tools::GFF->new(-file=>">$gfffile", -gff_version => 3);
  my @features = $self->db->features;
  foreach my $feat ( sort { $a->seq_id ge $b->seq_id } @features ){
    $gff->write_feature($feat);
  }
  return \@features;
}
sub write_genbank{
  my ($self,$outfile)=@_;
  my $gb = Bio::SeqIO->new(-file=>">$outfile", -format => 'genbank');
  my @seqids = $self->db->seq_ids;
  foreach my $seqid (@seqids){
    my $sequence = $self->db->fetch_sequence(-seq_id=>$seqid);
    if(1 > length $sequence){next;}
    my $seq = Bio::Seq::RichSeq->new(-id=>$seqid,-seq=>$sequence);
    my @feats = $self->db->get_features_by_location($seqid);
    $self->logmsg(sprintf("Writing %d features in %s\n",scalar @feats,$seqid));
    $seq->add_SeqFeature(@feats);
    $gb->write_seq($seq);
  }
}
sub load_gff{
  my ($self,$gfffile)=@_;
  my $gff = Bio::Tools::GFF->new(-file=>"<$gfffile", -gff_version => 3);
  my $count = 0;
  $self->logmsg("Loading $gfffile\n");
  while (my $feat = $gff->next_feature){
printf "Try to add feature %d:%d-%d\n",$count,$feat->start,$feat->end;
    $count++ if $self->add_feature($feat);
    if($count%1000 == 0){$self->logmsg("$count features loaded\r");}
  }
  $self->logmsg("$count features loaded\n");
}

sub load_fasta{
  my($self,$infile)=@_;
  my $fasta=Bio::SeqIO->new(-file=>"<$infile",-format=>'fasta')or die $!;
  $self->logmsg("Loading $infile\n");
  my $count = 0;
  while (my $seq = $fasta->next_seq){
   #$self->db->insert_sequence($seq->primary_id,$seq->seq,0);
    next unless $self->add_seq({ seq_id=>$seq->primary_id, seq=>$seq->seq});
    if($count%50 == 0){$self->logmsg("$count sequences loaded\r");}
    $count++;
  }
  $self->logmsg("$count sequences loaded\n");
}
1;
