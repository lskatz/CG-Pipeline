#!/usr/bin/env perl
use CGI qw/:standard/;
use DBI;
print header;
my $dbh=DBI->connect("DBI:mysql:host=localhost:3306;database=compgenomics","gbrowse2","Brows3r2") || die "cannot connect.$!";
my $sth=$dbh->prepare("select strain_name from genomeinfo_strain");
$sth->execute();
while($row=$sth->fetchrow_hashref()){
print "<td align=\"center\" >";
print "<b>".$row->{strain_name}."</b>";
push(@strains,$row->{strain_name});
print "</td>";
}
for($i=1;$i<=scalar(@strains);$i++){
$sth1=$dbh->prepare("select family from genomeinfo_moleculartyping where st_id in(select st_id from genomeinfo_strain_to_typing where strain_id=?) AND name='serogroup' ");
$sth1->execute($i);
while($row=$sth1->fetchrow_hashref()){
print $row->{family};
}
}
for($i=1;$i<=scalar(@strains);$i++){
$sth2=$dbh->prepare("select profile from genomeinfo_moleculartyping where st_id in (select st_id from genomeinfo_strain_to_typing where strain_id=?) AND profile IS NOT NULL");
$sth2->execute($i);
while($row1=$sth2->fetchrow_hashref()){
print $row1->{profile};
}
}
for($i=1;$i<=scalar(@strains);$i++){
$sth3=$dbh->prepare("select family from genomeinfo_moleculartyping where st_id in (select st_id from genomeinfo_strain_to_typing where strain_id=?) AND profile IS NOT NULL");
$sth3->execute($i);
while($row2=$sth3->fetchrow_hashref()){
print "ST-".$row2->{family};
}
}
$sth4=$dbh->prepare("select frequency_carriers from genomeinfo_strain");
$sth4->execute();
while($row3=$sth4->fetchrow_hashref()){
print $row3->{frequency_carriers};
}
$sth5=$dbh->prepare("select frequency_cases from genomeinfo_strain");
$sth5->execute();
while($row4=$sth5->fetchrow_hashref()){
print $row4->{frequency_cases};
}
$sth6=$dbh->prepare("select no_of_contigs from genomeinfo_strain");
$sth6->execute();
while($row5=$sth6->fetchrow_hashref()){
print $row5->{no_of_contigs};
}
$sth7=$dbh->prepare("select genome_size_bp from genomeinfo_strain");
$sth7->execute();
while($row6=$sth7->fetchrow_hashref()){
print $row6->{genome_size_bp};
}
$sth8=$dbh->prepare("select gc_content from genomeinfo_strain");
$sth8->execute();
while($row7=$sth8->fetchrow_hashref()){
print $row7->{gc_content};
}
$sth9=$dbh->prepare("select genbank_accno from genomeinfo_strain");
$sth9->execute();
while($row8=$sth9->fetchrow_hashref()){
print $row8->{genbank_accno};
}
$sth10=$dbh->prepare("select drs3_repeats from genomeinfo_strain");
$sth10->execute();
while($row9=$sth10->fetchrow_hashref()){
print $row9->{drs3_repeats};
}
$sth11=$dbh->prepare("select drs3_nf1_repeats from genomeinfo_strain");
$sth11->execute();
while($row10=$sth11->fetchrow_hashref()){
print $row10->{drs3_nf1_repeats};
}
for($i=1;$i<12;$i++){
$sth."$i"->disconnect();
}
