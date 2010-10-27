#!/usr/bin/env perl
package CGPBase;
use Exporter qw(import);
@EXPORT=qw(create_db sqlite2fasta makeblastdb);
use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;
sub create_db($$){
	my ($sqlitedb,$gfffile)=@_;
	my $cmd="bp_seqfeature_load.pl -a DBI::SQLite -c -d $sqlitedb $gfffile";
print STDERR "$cmd\n";
	system($cmd);
}
sub sqlite_connect($){
	my ($sqlitedb)=@_;
	return Bio::DB::SeqFeature::Store->new( -adaptor=>'DBI::SQLite',-dsn=>"$sqlitedb") or die "ERROR in CGPBase:$!\n";
}
sub sqlite2fasta{
	my ($sqlitedb,$outfile)=@_;
	my $db = sqlite_connect($sqlitedb);
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
sub makeblastdb($$$){
#args: organism id, fasta file, output dir
	my($organism,$file,$outdir)=@_;
    my $cmd="formatdb -t $organism -i $file -p F -n $outdir/$organism";
    system($cmd);
    $cmd="formatdb -t $organism -i $file -p T -n $outdir/$organism";
    system($cmd);
}
