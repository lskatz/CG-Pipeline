DROP TABLE IF EXISTS `blast_results`;
CREATE TABLE `blast_results` (
  `protein_id` varchar(50) NOT NULL,
  `uniprot_id` varchar(50) NOT NULL,
  `length` int(11) NOT NULL,
  `name` varchar(1000) default NULL,
  `rank` int(11) NOT NULL,
  `score` int(20) NOT NULL,
  `bits` float NOT NULL,
  `expectation` varchar(50) NOT NULL,
  `identity` float NOT NULL,
  `positives` float NOT NULL,
  KEY `protein_id` (`protein_id`),
  KEY `uniprot_id` (`uniprot_id`),
  KEY `rank` (`rank`)
);
DROP TABLE IF EXISTS `signalp_nn`;
CREATE TABLE `signalp_nn` (
  `protein_id` varchar(40) default NULL,
  `measure_type` varchar(16) default NULL,
  `position_start` int(20) default NULL,
  `position_end` int(20) default NULL,
  `value` double default NULL,
  `cutoff` double default NULL,
  `is_signal_peptide` varchar(5) default NULL,
  KEY (`protein_id`),
  KEY (`measure_type`)
);
DROP TABLE IF EXISTS `signalp_hmm`;
CREATE TABLE `signalp_hmm` (
  `protein_id` varchar(40) NOT NULL,
  `prediction` varchar(40) default NULL,
  `signal_peptide_probability` float default NULL,
  `max_cleavage_site_probability` float default NULL,
  `cleavage_start` int(10) default NULL,
  `cleavage_stop` int(10) default NULL,
  PRIMARY KEY (`protein_id`)
);
DROP TABLE IF EXISTS `tmhmm`;
CREATE TABLE `tmhmm` (
  `protein_id` varchar(40) NOT NULL default '',
  `predicted_number` int(10) default NULL,
  `length` int(11) default NULL,
  `expected_number_aa` double default NULL,
  `expected_number_aa_60` double default NULL,
  `total_prob_n_in` double default NULL,
  PRIMARY KEY  (`protein_id`)
);
DROP TABLE IF EXISTS `tmhmm_location`;
CREATE TABLE `tmhmm_location` (
  `protein_id` varchar(50) NOT NULL,
  `location` varchar(50) NOT NULL,
  `from_pos` int(20) default NULL,
  `to_pos` int(20) default NULL
);
DROP TABLE IF EXISTS `uniprot`;
CREATE TABLE `uniprot` (
  `accession` varchar(50) NOT NULL,
  `name` varchar(500) default NULL,
  `dataset` varchar(50) NOT NULL,
  `protein_name` varchar(1000) default NULL,
  `protein_type` varchar(1000) default NULL,
  `gene_type` varchar(500) default NULL,
  `gene_name` varchar(1000) default NULL,
  `gene_id` int(11) default NULL,
  PRIMARY KEY  (`accession`)
);
DROP TABLE IF EXISTS `uniprot_evidence`;
CREATE TABLE `uniprot_evidence` (
  `accession` varchar(50) NOT NULL,
  `evidence_id` varchar(50) NOT NULL,
  `database_name` varchar(500) NOT NULL,
  `name` varchar(5000) default NULL,
  KEY `accession` (`accession`),
  KEY `evidence_id` (`evidence_id`)
);
DROP TABLE IF EXISTS `ipr_hits`;
CREATE TABLE `interpro_hits` (
  `protein_id` varchar(50) NOT NULL,
  `length` integer,
  `domain_id` varchar(50) NOT NULL default 'noIPR',
  `name` varchar(1000) default NULL,
  `type` varchar(100) default NULL,
   KEY  (`protein_id`),
   KEY  (`domain_id`)
);
DROP TABLE IF EXISTS `ipr_matches`;
CREATE TABLE `ipr_matches` (
  `protein_id` varchar(50) NOT NULL,
  `accession_num` varchar(50) default NULL,
  `database_name` varchar(100) NOT NULL,
  `start` int(20) default NULL,
  `end` int(20) default NULL,
  `evalue` varchar(40) default NULL,
  `status` varchar(5) NOT NULL,
  `evidence` varchar(100) NOT NULL,
  KEY `protein_id` (`protein_id`),
  KEY `accession_num` (`accession_num`),
  KEY `evidence` (`evidence`)
);
DROP TABLE IF EXISTS `ipr_classifications`;
CREATE TABLE `ipr_classifications` (
  `protein_id` varchar(50) NOT NULL,
  `go_id` varchar(50) NOT NULL default '',
  `class_type` varchar(16) default NULL,
  `category` varchar(50) default NULL,
  `description` varchar(2000) default NULL,
  KEY  (`protein_id`),
  KEY  (`go_id`)
);
DROP TABLE IF EXISTS `ipr_child_references`;
CREATE TABLE `ipr_child_references` (
  `protein_id` varchar(50) NOT NULL,
  `interpro_hit_id` varchar(50) NOT NULL,
  `interpro_child_reference` varchar(50) NOT NULL
);
DROP TABLE IF EXISTS `vfdb_hits`;
CREATE TABLE `vfdb_hits` (
  `protein_id` varchar(50) NOT NULL,
  `target_id` varchar(50) NOT NULL,
  `expectation` float,
  `coverage` float,
  `db_name` varchar(50) default NULL,
  `identity` float,
  KEY (`protein_id`),
  KEY (`target_id`)
);
