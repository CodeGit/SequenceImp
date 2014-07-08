#!/usr/bin/perl

use warnings;
use strict;
use Test::More;
use File::Spec;
use File::Basename;
use File::Path;
use File::Copy;
use Env::Path;
use IPC::Open3;
use Cwd 'abs_path';
use Data::Dumper;

my $SEQIMP = "imp_commandline.pl";
my $CWD = dirname(abs_path($0));
my $SEQIMP_ROOT = abs_path(File::Spec->catfile($CWD, ".."));
my $DEPENDENCIES_PATH = File::Spec->catfile($SEQIMP_ROOT, "dependencies-bin");
my $TEST_DIR = File::Spec->catfile($SEQIMP_ROOT, "test-data");
my $TEST_SEQUENCE_STAGING = File::Spec->catfile($TEST_DIR, "seq-staging");
my $TEST_METADATA = File::Spec->catfile($TEST_DIR, "metadata");
my $TEST_RUN = File::Spec->catfile($TEST_DIR, "run");
my $TEST_ANALYSIS = File::Spec->catfile($TEST_RUN, "analysis");
my $ANNOTATION = File::Spec->catfile($SEQIMP_ROOT, "Annotation-12-164");

my @DESCRIPTION_COLUMNS = qw(Name File Geometry Barcodes 5p_ad 3p_ad 5p_seq_insert 3p_seq_insert);

my @CONFIG_ENTRIES;
push(@CONFIG_ENTRIES, qw(@organise reap2imp geoConfig));
push(@CONFIG_ENTRIES, qw(@reaper+filter fastq plotZeros));
push(@CONFIG_ENTRIES, qw(@reaper reapConfig perlUnq));
push(@CONFIG_ENTRIES, qw(@filter low minSize maxSize five three));
push(@CONFIG_ENTRIES, qw(@align+features genome ensversion));
push(@CONFIG_ENTRIES, qw(@align chunk mismatches maxHits sam));
push(@CONFIG_ENTRIES, qw(@features feature mirversion annot_conflict overlap proportional separate_loci repMaxHits repMismatches repChunk));

$ENV{'SEQIMP_ROOT'} = $SEQIMP_ROOT;
my $path = Env::Path->PATH;
$path->Prepend(File::Spec->catfile($SEQIMP_ROOT, "bin"));
$path->Prepend($DEPENDENCIES_PATH);
$path->Prepend(File::Spec->catfile($DEPENDENCIES_PATH, "R", "bin"));


sub cleanup {
	my $dir = @_;
	rmtree($TEST_DIR);
}

sub prepare {
	mkpath($TEST_DIR);
	mkpath($TEST_SEQUENCE_STAGING);
	mkpath($TEST_METADATA);
	mkpath($TEST_RUN);
}

sub writeDescription {
	my ($description, $path) = @_;
	my $outFile = File::Spec->catfile($TEST_METADATA, "description.txt");
	
	open(OUT, ">$outFile") or die ("Failed to open $outFile for writing:$!");
	
	print OUT join("\t", @DESCRIPTION_COLUMNS);
	print OUT "\n";
	for my $line (@{$description}) {
		for my $col (@DESCRIPTION_COLUMNS) {
			my $value = $line->{$col};
			$value = "-" unless (defined $value);
			print OUT "$value";
			if ($col eq $DESCRIPTION_COLUMNS[$#DESCRIPTION_COLUMNS]) {
				print OUT "\n";
			} else {
				print OUT "\t";	
			} 
		}
	}
	close(OUT);
	return $outFile;
}

sub writeConfig {
	my ($config, $path) = @_;
	
	my $outFile = File::Spec->catfile($TEST_METADATA, "config.txt");
	
	open(OUT, ">$outFile") or die ("Failed to open $outFile for writing:$!");
	foreach my $line (@CONFIG_ENTRIES) {
		if ($line =~ /^\@/) {
			print OUT "$line\n";
		} else {
			my $value = $config->{$line};
			if (defined $value) {
				print OUT "$line\t$value\n";
			} else {
				print OUT "$line\tNA\n";
			}
		}
	}
	close(OUT);
	return $outFile;
}

sub runOrganise {
	my ($dataDir, $analysisDir, $configFile, $descriptionFile) = @_;
	
	my $pid = open3(\*IN, \*OUT, \*ERR, 
	"$SEQIMP --step=organise --no-unique --description=$descriptionFile --user-configuration=$configFile --dataDir=$dataDir --outDir=$analysisDir");
	close(IN);
	my @outlines = <OUT>;
	my @errlines = <ERR>;
	waitpid($pid, 0);
	
	is ($?, 0, "Organise exit code");
	print "ERROR\n@errlines\n" if ($?);	
}

sub runReaper{
	my ($dataDir, $analysisDir, $configFile, $descriptionFile) = @_;
	
	$analysisDir = File::Spec->catfile($analysisDir, "analysis");
	my $pid = open3(\*IN, \*OUT, \*ERR, 
	"$SEQIMP --step=reaper --description=$descriptionFile --user-configuration=$configFile --dataDir=$dataDir --analysisDir=$analysisDir --annotationDir=$ANNOTATION");
	close(IN);
	my @outlines = <OUT>;
	my @errlines = <ERR>;
	waitpid($pid, 0);
	
	is ($?, 0, "Filter exit code");
	print "ERROR\n@errlines\n" if ($?);
}

sub runSeqImp {
	my ($dataDir, $analysisDir, $configFile, $descriptionFile) = @_;
	
	runOrganise($dataDir, $analysisDir, $configFile, $descriptionFile);
    runReaper($dataDir, $analysisDir, $configFile, $descriptionFile);
    #runFilter($dataDir, $analysisDir, $configFile, $descriptionFile);
    #runAlign($dataDir, $analysisDir, $configFile, $descriptionFile);
    #runFeatures($dataDir, $analysisDir, $configFile, $descriptionFile);
}

sub test_initialise {
	my $pid = open3(\*IN, \*OUT, \*ERR, "$SEQIMP --system-check");
	close(IN);
	my @outlines = <OUT>;
	my @errlines = <ERR>;
	waitpid($pid, 0);
	
	is ($?, 0, "system check exit code"); 
}

sub test_five_prime_barcode_run {
	cleanup();
	prepare();
	
	my $seqfile = "testfile.fastq.gz";
	my $seqLocation = File::Spec->catfile($CWD, "resources", "5p_barcode_and_insert_data");
	my $seqfileLocation = File::Spec->catfile($seqLocation, $seqfile);
	my $newLocation = File::Spec->catfile($TEST_SEQUENCE_STAGING, $seqfile);
	copy($seqfileLocation, $newLocation);
	#ok(-e $seqfileLocation, "Test sequence file");
	#ok(-e $newLocation, "Test sequence file copy to '$newLocation'");
	
	my $description = [];
	my $file = {};
	#->{'analysis'} = "test5p";
	$file->{"Name"} = "testfile";
	$file->{"File"} = $seqfile;
	$file->{'Geometry'} = "5p_barcode_and_insert";
	$file->{'Barcodes'} = "CATG,CTAG,GATC,ACTG,GCTA,GTAC,AGTC,GTCA,CTGA,TCAG";
	$file->{'5p_ad'} = "GTTCAGAGTTCTACAGTCCGACGATC";
	$file->{'3p_ad'} = "ATCTCGTATGCCGTCTTCTGCTTG";
	$file->{'5p_seq_insert'} = "CG";
	$file->{'3p_seq_insert'} = "-";
	push(@{$description}, $file);
	my $descriptionFile = writeDescription($description, $TEST_SEQUENCE_STAGING);
	
	my $config = {};
	$config->{'minSize'} = 18;
	$config->{'maxSize'} = 26;
	$config->{'genome'} = "mouse";
	$config->{'ensversion'} = 66;
	$config->{'mismatches'} = 2;
	$config->{'maxHits'} = 20;
	$config->{'sam'} = "FLAG";
	$config->{'feature'} = "miRNA";
	$config->{'mirversion'} = 18;
	$config->{'annot_conflict'} = "merge";
	$config->{'overlap'} = 15;
	my $configFile = writeConfig($config, $TEST_SEQUENCE_STAGING);
	
    runSeqImp($TEST_SEQUENCE_STAGING, $TEST_RUN, $configFile, $descriptionFile);
    
    #organise tests
    ok(-e $TEST_ANALYSIS, "Organise analysis directory created");
    my $testfileAnalysis = File::Spec->catfile($TEST_ANALYSIS, "testfile");
    ok(-e $testfileAnalysis, "Organise testfile subdirectory created");
    my $testfileMetadata = File::Spec->catfile($testfileAnalysis, "metadata");
    ok(-e $testfileMetadata, "Organise testfile metadata subdirectory created");
    my $testMetadataFile = File::Spec->catfile($testfileMetadata, "metadata.txt");
    ok(-e $testfileMetadata, "Organise testfile metadata file created");
    my $testfileData = File::Spec->catfile($testfileAnalysis, "data");
    ok(-e $testfileData, "Organise testfile data subdirectory created");
    my $testfileDataFile = File::Spec->catfile($testfileData, "analysis_$seqfile");
    ok(-e $testfileDataFile, "Organise testfile data file created");
    #reaper tests
    my $testfileReaper = File::Spec->catfile($testfileAnalysis, "REAPER");
    ok(-e $testfileReaper, "Reaper analysis directory created");
    #barcode files
	for my $barcode (split(",", $file->{'Barcodes'})) {
    	for my $suffix (qw(clean.gz report.input.q report.input.nt report.clean.nt report.clean.len uniquify.log report.clean.trinucl report.clean.annotlen clean.uniquified.fa.gz)) {
    		my $testReaperBarcodeOutputFile = File::Spec->catfile($testfileReaper, "testfile.$barcode.$suffix");
    		ok(-e $testReaperBarcodeOutputFile, "Reaper analysis file testfile.$barcode.$suffix created");
    	}
    }
    #summary files
    for my $suffix (qw(total.report.miss.nt total.report.miss.len total.report.input.q total.report.input.nt sumstat lint.gz)) {
    	my $testReaperTotalOutputFile = File::Spec->catfile($testfileReaper, "testfile.$suffix");
    	ok(-e $testReaperTotalOutputFile, "Reaper analysis file testfile.$suffix created");
    }
    my $testfileQC = File::Spec->catfile($testfileAnalysis, "QC");
    ok(-e $testfileQC, "Reaper analysis QC directory created");
    my $testfileQCFile = File::Spec->catfile($testfileQC, "testfile_Reaper_qc.pdf");
    ok(-e $testfileQCFile, "Reaper analysis QC pdf file created");
    
}

test_initialise();
test_five_prime_barcode_run();

done_testing();
