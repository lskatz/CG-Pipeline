# $Id: $

=head1 NAME

CGPipeline::Store::Prediction - plugin for Prediction pipeline functions 

=head1 AUTHOR

Jay Humphrey <jhumphrey6@gmail.com>

=cut

=head1 SYNOPSIS

=head1 TODO

finish perldoc...

=cut

package CGPipeline::Store::Prediction; 
use strict;
use String::Scanf;
use base qw/CGPipeline::Store/;

=cut

head2  new
 Title   : new
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    : 

=cut

sub new{
  my $self=shift;
  my %params = @_;
  my $object = $self->SUPER::new(@_);
  my $locidfmt = $params{locus_id_format};
  $object->{locus_id_format} = $locidfmt || "LOC%6d";
  return $object;
}

=cut

head2  clean
 Title   : clean
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    : 

=cut

sub clean{
  my $self=shift;
  $self->db->delete($self->db->features);
}
=cut

head2  add_feature
 Title   : add_feature
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    : 

=cut
sub add_feature{ 
  my $self=shift; #feature should be a Bio::SeqFeatureI
  my $newfeat=shift;
  my ($seq_id,$locus_id) = ($newfeat->seq_id,"");
  ($locus_id) = $newfeat->get_tag_values('locus_id') if $newfeat->has_tag('locus_id');
  if(!$locus_id){
    $locus_id = $self->next_locus_id;
    $newfeat->add_tag_value('locus_id',$locus_id);
  }
  if(!$self->seq_exists($seq_id)){
    $self->logmsg("Warning:sequence $seq_id does not exist\n");
  }
  my @features = $self->db->features(-attributes=>{'locus_id'=>$locus_id});#,-seq_id=>$seq_id);
  push(@features,$self->db->get_features_by_location($seq_id,$newfeat->start,$newfeat->end));
  foreach my $feature(@features){
    if($self->{volatile}){
      $self->logmsg("Replacing feature $seq_id:$locus_id\n");
      $self->db->delete(($feature));
    }
    else{
      $self->logmsg("\t\t\t\tskipping $locus_id (exists)\n");
      return 0;
    }
  }
  $self->db->store($newfeat) or die "Error:failed at storing feature:$!\n";
  
  return 1;
}

=cut

head2  next_locus_id
 Title   : next_locus_id
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    : 

=cut
sub next_locus_id{
  my $self=shift;
  my $format = $self->{'locus_id_format'};
  my $id=$self->{'lastID'} || 1;
  my $locus_id = "";
  $self->db->index_tables;
  $self->db->optimize;
  while($locus_id eq ""){
    my $newid = sprintf($format,$id);
    my ($feature)= $self->db->features(-attributes=>{'locus_id'=>$newid});
    if($feature){$id++;}
    else{$locus_id = $newid;last;} 
  }
  $self->{'lastID'}=$id;
  return $locus_id;
}
1;
