# $Id: $

=head1 NAME

CGPipeline::Store::Export - plugin for Export pipeline functions 

=head1 AUTHOR

Jay Humphrey <jhumphrey6@gmail.com>

=cut

=head1 SYNOPSIS

=head1 TODO

finish perldoc...

=cut

package CGPipeline::Store::Export; 
use base qw/CGPipeline::Store/;
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
    my $stable_id = ($ftr->get_tag_values('stable_id'))[0];
    my $defline = sprintf("lcl|%s|%s|%d|%d",$stable_id,$ftr->primary_tag,$ftr->start,$ftr->end);
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
    my $stable_id = ($ftr->get_tag_values('stable_id'))[0];
    my $defline = sprintf("lcl|%s|%s|%d|%d",$stable_id,$ftr->primary_tag,$ftr->start,$ftr->end);
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
