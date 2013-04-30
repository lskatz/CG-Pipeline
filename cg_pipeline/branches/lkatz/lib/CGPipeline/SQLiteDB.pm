#!/usr/bin/env perl
package CGPipeline::SQLiteDB; 
use Bio::SeqIO;
use Data::Dumper;
use Bio::DB::SeqFeature::Store;
use Bio::Tools::GFF;
sub new{
	my $class = shift;
	my %params = @_;
	die "sqlite database file or new filename must be defined\n" if !defined $params{'file'};
	$params{'verbose'}=0 unless defined $params{'verbose'};
	my $object = bless({
		file=>$params{'file'},
		verbose=>$params{'verbose'},
		db=>undef
	},$class);
	$object->_connect;
	return $object;
}
sub logmsg{
	my $self = shift;
#if(! $self->{'verbose'} ){return;}
	my $msg = shift;
	printf STDERR "$msg";
}
sub db{
	my $self = shift;
	return $self->{'db'};
}
sub _connect($){
	my $self=shift;
	my $is_new = 0;
	my $sqlitedb = $self->{'file'};
	if( ! -e $sqlitedb){$is_new=1;}
	$self->{'db'} = Bio::DB::SeqFeature::Store->new( -adaptor=>'DBI::SQLite',-dsn=>"$sqlitedb", -create=>$is_new) or die "ERROR in CGPBase:$!\n";
}
sub write_gff{
	my ($self,$gfffile)=@_;
	my $db = $self->db;
	my $gff = Bio::Tools::GFF->new(-file=>">$gfffile", -gff_version => 3);
	my @features = $db->features;
	foreach my $feat ( sort { $a->seq_id ge $b->seq_id } @features ){
		$gff->write_feature($feat);
	}
}
sub write_genbank{
	my ($self,$outfile)=@_;
	my $db = $self->db;
	my $gb = Bio::SeqIO->new(-file=>">$outfile", -format => 'genbank');
	my @seqids = $db->seq_ids;
	foreach my $seqid (@seqids){
		my $sequence = $db->fetch_sequence(-seq_id=>$seqid);
		if(1 > length $sequence){next;}
		my $seq = Bio::Seq::RichSeq->new(-id=>$seqid,-seq=>$sequence);
		my @feats = $db->get_features_by_location($seqid);
		$self->logmsg(sprintf("Writing %d features in %s\n",scalar @feats,$seqid));
		$seq->add_SeqFeature(@feats);
		$gb->write_seq($seq);
	}
}
sub load_gff{
	my ($self,$gfffile)=@_;
	my $db = $self->db;
	my $gff = Bio::Tools::GFF->new(-file=>"<$gfffile", -gff_version => 3);
	my $count = 0;
	$self->logmsg("Loading $gfffile\n");
	while (my $feat = $gff->next_feature){
		$db->store($feat);
		if($count%1000 == 0){$self->logmsg("$count features loaded\r");}
		$count++;
	}
	$self->logmsg("$count features loaded\n");
}

sub load_fasta{
	my($self,$infile)=@_;
	my $db = $self->db;
	my $fasta=Bio::SeqIO->new(-file=>"<$infile",-format=>'fasta')or die $!;
	$self->logmsg("Loading $infile\n");
	my $count = 0;
	while (my $seq = $fasta->next_seq){
		$db->insert_sequence($seq->primary_id,$seq->seq,0);
		if($count%1000 == 0){$self->logmsg("$count sequences loaded\r");}
		$count++;
	}
	$self->logmsg("$count sequences loaded\n");
}
sub get_contigs{
	my ($self,$outfile)=@_;
	my $db = $self->db;
	my $fasta=Bio::SeqIO->new(-file=>">$outfile",-format=>'fasta')or die $!;
	my $seqio = $db->get_seq_stream;
	while (my $seq = $seqio->next_seq ){
		$fasta->write_seq($seq->seq);
	}
}
sub get_fasta{
	my ($self,$outfile)=@_;
	my $db = $self->db;
	my $fasta=Bio::SeqIO->new(-file=>">$outfile",-format=>'fasta')or die $!;
	my @features=$db->features(-type=>'CDS');
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
	my $db = $self->db;
	my $fasta=Bio::SeqIO->new(-file=>">$outfile",-format=>'fasta')or die $!;
	my @features=$db->features(-type=>'CDS');
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
sub makeblastdb($$$$){
#args: organism id, fasta file, protein(T/F), output dir
	my($organism,$file,$prot,$outdir)=@_;
    my $cmd="formatdb -t $organism -i $file -p $prot -n $outdir/$organism";
    system($cmd);
}
1;
