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
use Bio::SeqIO;
use Data::Dumper;
use Bio::DB::SeqFeature::Store;
use Bio::Tools::GFF;
sub new{
  my $class = shift;
  my %params = @_;
  die "new or existing sqlite database file must be defined\n" unless defined $params{'file'};
  $paramsverbose=0 unless defined $paramsverbose;
  my $object = bless({
    file=>$params{'file'},
    verbose=>$paramsverbose,
    volatile=>$paramsvolatile || 0,
    db=>undef
  },$class);
  $object->_connect;
  return $object;
}
sub logmsg{
  my $self = shift;
#if(! $self->verbose ){return;}
  my $msg = shift;
  printf STDERR "$msg";
}
sub db{
  my $self = shift;
  return $self->{'db'};
}
sub verbose{
  my $self = shift;
  return $self->verbose;
}
sub volatile{
  my $self = shift;
  return $self->volatile;
}
sub _connect{
  my $self=shift;
  my $is_new = 0;
  my $sqlitedb = $self->{'file'};
  if( ! -e $sqlitedb){$is_new=1;}
  $self->{'db'} = Bio::DB::SeqFeature::Store->new( -adaptor=>'DBI::SQLite',-dsn=>"$sqlitedb", -create=>$is_new);
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
  if( my $stored_seq = $self->get_seq_by_id($seq->{'seq_id'})){
    logmsg(sprintf("Warning:%s already exists\n",$seq->{'seq_id'}));
    return $stored_seq unless $self->volatile;
  }
  if( ! $self->sane_id($seq->{'seq_id'}){
    logmsg(sprintf("Failed:%s is unsafe\n",$seq->{'seq_id'}));
    return 0;
  }
  return $self->db->insert_sequence($seq->{'seq_id'},$seq->{'seq'},0);
}
sub get_seq_by_id{
  my ($self,$seqid)=@_;
  return $self->db->fetch_sequence(-seq_id=>$seqid,-bioseq=>1);
}
sub add_locus{ #belongs in CGPipeline::SQLiteDB::Predict ?
  my ($self,$params)=@_; #feature should be a Bio::SeqFeatureI
  my $seq_id = $params->{seq_id};
  my $name = $params->{name};
  my $seq=$self->get_seq_by_id($seq_id, -bioseq=>1);
  if(!$seq){
    logmsg("Warning:sequence $seq_id does not exist\n");
  }
  my ($feature) = $self->db->features(-name=>$name,-seq_id=>$seq_id);
  if($feature){
    if($self->volatile){
      logmsg("Replacing feature $seq_id:$name\n");
      $self->db->delete(($feature));
    }
    else{return $feature;}
  }
  my $feature = Bio::SeqFeature::Generic->new(@$params);
  return $self->db->store($feature);
}
sub add_feature{ # all features will be children of loci
  
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
    $self->db->store($feat);
    if($count%1000 == 0){$self->logmsg("$count features loaded\r");}
    $count++;
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
    next unless $self->add_seq(-seq_id=>$seq->primary_id,-seq=>$seq->seq);
    if($count%1000 == 0){$self->logmsg("$count sequences loaded\r");}
    $count++;
  }
  $self->logmsg("$count sequences loaded\n");
}
sub get_contigs{
  my ($self,$outfile)=@_;
  my $fasta=Bio::SeqIO->new(-file=>">$outfile",-format=>'fasta')or die $!;
  my $seqio = $self->db->get_seq_stream;
  while (my $seq = $seqio->next_seq ){
    $fasta->write_seq($seq->seq);
  }
}
sub get_fasta{
  my ($self,$outfile)=@_;
  my $fasta=Bio::SeqIO->new(-file=>">$outfile",-format=>'fasta')or die $!;
  my @features=$self->db->features(-type=>'CDS');
  foreach my $ftr(@features){
    my $seq2fasta=$ftr->seq;
    my $locus_tag = ($ftr->get_tag_values('locus_tag'))[0];
    my $defline = sprintf("lcl|%s|%s|%d|%d",$locus_tag,$ftr->primary_tag,$ftr->start,$ftr->end);
    $seq2fasta->display_name($defline);
    if($ftr->has_tag('gene')){$seq2fasta->desc(($ftr->get_tag_values('gene'))[0]);}
    elsif($ftr->has_tag('product')){$seq2fasta->desc(($ftr->get_tag_values('product'))[0]);}
    elsif($ftr->primary_tag =~ /gene/){}#skip it! avoids printing duplicate gene/CDS pairs.
    else{$seq2fasta->desc("predicted cds");}
    my $seq=$seq2fasta->seq;
    $seq=~s/-/N/g;
    $seq2fasta->seq($seq);
    $fasta->write_seq($seq2fasta) or die "$!\n";
  }
  return $outfile;
}
sub get_fastaprot{
  my ($self,$outfile)=@_;
  my $fasta=Bio::SeqIO->new(-file=>">$outfile",-format=>'fasta')or die $!;
  my @features=$self->db->features(-type=>'CDS');
  foreach my $ftr(@features){
    my $seq2fasta=$ftr->seq;
    my $locus_tag = ($ftr->get_tag_values('locus_tag'))[0];
    my $defline = sprintf("lcl|%s|%s|%d|%d",$locus_tag,$ftr->primary_tag,$ftr->start,$ftr->end);
    $seq2fasta->display_name($defline);
    if($ftr->has_tag('gene')){$seq2fasta->desc(($ftr->get_tag_values('gene'))[0]);}
    elsif($ftr->has_tag('product')){$seq2fasta->desc(($ftr->get_tag_values('product'))[0]);}
    elsif($ftr->primary_tag =~ /gene/){}#skip it! avoids printing duplicate gene/CDS pairs.
    else{$seq2fasta->desc("predicted cds");}
    my $seq=$seq2fasta->translate->seq;
    $seq=~s/-/N/g;
    $seq2fasta->seq($seq);
    $fasta->write_seq($seq2fasta) or die "$!\n";
  }
  return $outfile;
}
1;
