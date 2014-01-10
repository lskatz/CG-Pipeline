#!/usr/bin/env perl
# TODO multithread SeqIO/in
# TODO buffer the output
# TODO download the correct genus files
# wget http://www.bacterio.net/-ac.html http://www.bacterio.net/-dl.html http://www.bacterio.net/-mr.html http://www.bacterio.net/-sz.html
# mv ./-ac.html ac.html
# mv ./-dl.html dl.html
# mv ./-mr.html mr.html
# mv ./-sz.html sz.html


use Bio::Perl;
use Data::Dumper;

sub logmsg{print STDERR "@_\n";}
exit main();

sub main{
  open(GENUS,"genus.txt");
  while(<GENUS>){
    chomp; 
    $genus{$_}=1; 
    push(@genus,$_); 
  } 
  close GENUS; 
  $out=Bio::SeqIO->new(-format=>"fasta"); 
  $in=Bio::SeqIO->new(-file=>"uniprot_sprot_trembl.fasta"); 
  my $genusRegex=join("|",@genus);
  while($seq=$in->next_seq){
    logmsg $i if(++$i % 1000 == 0); 
    $id=$seq->id." ".$seq->desc; 
    if($id=~/$genusRegex/){
      $out->write_seq($seq);
    }
  }
  
  return 0;
}
