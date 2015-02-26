#!/usr/bin/env perl

###########################################################################################################
# Integration test of seqimp covers all steps from organise to features
###########################################################################################################

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
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PathUtil qw(findExecutable);

my $SEQIMP = "imp_commandline.pl";

my $CWD = dirname(abs_path($0));
my $SEQIMP_ROOT = abs_path(File::Spec->catfile($CWD, ".."));
my $DEPENDENCIES_PATH = File::Spec->catfile($SEQIMP_ROOT, "dependencies-bin");
my $OSPATH = getOSPath($DEPENDENCIES_PATH);
my $RPATH = File::Spec->catfile($OSPATH, "R", "bin");
my $PERLPATH = File::Spec->catfile($OSPATH, "perl", "bin");

my $TEST_DIR = File::Spec->catfile($SEQIMP_ROOT, "test-data");
my $TEST_SEQUENCE_STAGING = File::Spec->catfile($TEST_DIR, "seq-staging");
my $TEST_METADATA = File::Spec->catfile($TEST_DIR, "metadata");
my $TEST_RUN = File::Spec->catfile($TEST_DIR, "run");
my $TEST_ANALYSIS = File::Spec->catfile($TEST_RUN, "analysis");
my $ANNOTATION = File::Spec->catfile($SEQIMP_ROOT, "annotation_210814");

my @DESCRIPTION_COLUMNS = qw(Name File Geometry Barcodes 5p_ad 3p_ad 5p_seq_insert 3p_seq_insert);

my @CONFIG_ENTRIES;
push(@CONFIG_ENTRIES, qw(@organise reap2imp geoConfig));
push(@CONFIG_ENTRIES, qw(@reaper+filter fastq plotZeros));
push(@CONFIG_ENTRIES, qw(@reaper reapConfig perlUnq));
push(@CONFIG_ENTRIES, qw(@filter low minSize maxSize five three));
push(@CONFIG_ENTRIES, qw(@align+features genome ensversion));
push(@CONFIG_ENTRIES, qw(@align chunk mismatches maxHits sam));
push(@CONFIG_ENTRIES, qw(@features feature mirversion annot_conflict overlap proportional separate_loci collapse_method repversion repMaxHits repMismatches repChunk));

$ENV{'SEQIMP_ROOT'} = $SEQIMP_ROOT;
my $path = Env::Path->PATH;
$path->Prepend(File::Spec->catfile($SEQIMP_ROOT, "bin"));
$path->Prepend($OSPATH);
$path->Prepend($RPATH);
$path->Prepend($PERLPATH);

sub getOSPath {
	my ($dependenciesPath) = @_;
	my $osPath;
	my $os = $^O;
	if ($os eq "darwin") {
		$osPath = File::Spec->catfile($dependenciesPath, "macosx", "bin");
	} elsif ($os eq "MSWin32"){
		$osPath = File::Spec->catfile($dependenciesPath, "windows", "bin");
	} elsif ($os eq "linux") {
		$osPath = File::Spec->catfile($dependenciesPath, "linux", "bin");
	}
	return $osPath;
}

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

sub runStep {
	my ($step, $configFile, $descriptionFile, $args) = @_;

	my $pid = open3(\*IN, \*OUT, \*ERR, 
	"$SEQIMP --debug --step=$step --description=$descriptionFile --user-configuration=$configFile $args");
	close(IN);
	my @outlines = <OUT>;
	my @errlines = <ERR>;
	waitpid($pid, 0);
	
	is ($?, 0, "$step exit code");
	print "ERROR\n@errlines\n" if ($?);
}

sub test_Prerequisites {
	ok (-e $OSPATH, "Found dependencies");
	ok (-e $RPATH, "Found R bin directory");
	ok (-e $PERLPATH, "Found perl bin directory");
	ok (-e $ANNOTATION, "Found annotation directory");
	
	my @executables = qw(perl R swan tally samtools reaper minion bowtie);
	foreach my $executable (@executables) {
		my $executablePath = findExecutable($executable);
		ok (-e $executablePath, "Found $executable at $executablePath");
	}
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
	ok(-e $seqfileLocation, "Test sequence file");
	ok(-e $newLocation, "Test sequence file copy to '$newLocation'");
	my $analysisDir = File::Spec->catfile($TEST_RUN, "analysis");
	
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
	$config->{'ensversion'} = 1;
	$config->{'mismatches'} = 2;
	$config->{'maxHits'} = 20;
	$config->{'sam'} = "FLAG";
	$config->{'feature'} = "miRNA";
	$config->{'mirversion'} = 1;
	$config->{'annot_conflict'} = "merge";
	$config->{'overlap'} = 15;
	#newer options
	$config->{'proportional'} = "NA";
	$config->{'separate_loci'} = "NA";
	$config->{'collapse_method'} = "mature_id";
	$config->{'repversion'} = "NA";
	$config->{'repMaxHits'} = "NA";
	$config->{'repMismatches'} = "NA";
	$config->{'repChunk'} = "NA";
	
	my $configFile = writeConfig($config, $TEST_SEQUENCE_STAGING);
	
    #runSeqImp($TEST_SEQUENCE_STAGING, $TEST_RUN, $configFile, $descriptionFile);
    
    #organise tests
    runStep("organise", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --outDir=$TEST_RUN");
    ok(-e $TEST_ANALYSIS, "Organise directory created");
    my $testfileAnalysis = File::Spec->catfile($TEST_ANALYSIS, "testfile");
    ok(-e $testfileAnalysis, "Organise seqfile subdirectory created");
    my $testfileMetadata = File::Spec->catfile($testfileAnalysis, "metadata");
    ok(-e $testfileMetadata, "Organise seqfile metadata subdirectory created");
    my $testMetadataFile = File::Spec->catfile($testfileMetadata, "metadata.txt");
    ok(-e $testfileMetadata, "Organise seqfile metadata file created");
    my $testfileData = File::Spec->catfile($testfileAnalysis, "data");
    ok(-e $testfileData, "Organise seqfile data subdirectory created");
    my $testfileDataFile = File::Spec->catfile($testfileData, "analysis_$seqfile");
    ok(-e $testfileDataFile, "Organise seqfile data file created");
    
    #reaper tests
    runStep("reaper", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir");
    my $testfileReaper = File::Spec->catfile($testfileAnalysis, "REAPER");
    ok(-e $testfileReaper, "Reaper analysis directory created");
    #barcode files
	for my $barcode (split(",", $file->{'Barcodes'})) {
    	for my $suffix (qw(clean.gz report.input.q report.input.nt report.clean.nt report.clean.len uniquify.log report.clean.trinucl report.clean.annotlen clean.uniquified.fa.gz)) {
    		my $testReaperBarcodeOutputFile = File::Spec->catfile($testfileReaper, "testfile.$barcode.$suffix");
    		ok(-e $testReaperBarcodeOutputFile, "Reaper file testfile.$barcode.$suffix created");
    	}
    }
    #summary files
    for my $suffix (qw(total.report.miss.nt total.report.miss.len total.report.input.q total.report.input.nt sumstat lint.gz)) {
    	my $testReaperTotalOutputFile = File::Spec->catfile($testfileReaper, "testfile.$suffix");
    	ok(-e $testReaperTotalOutputFile, "Reaper file testfile.$suffix created");
    }
    my $testfileQC = File::Spec->catfile($testfileAnalysis, "QC");
    ok(-e $testfileQC, "Reaper analysis QC directory created");
    my $testfileReaperQCFile = File::Spec->catfile($testfileQC, "testfile_Reaper_qc.pdf");
    ok(-e $testfileReaperQCFile, "Reaper QC pdf file created");
    
    #filter tests
    runStep("filter", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir");
    my $testfileFilter = File::Spec->catfile($testfileAnalysis, "PROCESSED");
    ok(-e $testfileFilter, "Filter directory created");
    my $testfileFilterTrimmed = File::Spec->catfile($testfileFilter, "total.saved.trimmit.tab");
    ok(-e $testfileFilter, "Filter trimmed file created");
    my $testfileFilterQC = File::Spec->catfile($testfileQC, "testfile_Processed_reads_qc.pdf");
    ok(-e $testfileFilterQC, "Filter QC file created");
    #barcode files
	for my $barcode (split(",", $file->{'Barcodes'})) {
		for my $suffix (qw(report.clean.processed.annotlen clean.processed.fa.gz)) {
    		my $testFilterBarcodeOutputFile = File::Spec->catfile($testfileFilter, "testfile.$barcode.$suffix");
    		ok(-e $testFilterBarcodeOutputFile, "Filter file testfile.$barcode.$suffix created");
    	}
	}
	
	#align tests
	runStep("align", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir --annotationDir=$ANNOTATION");
    my $testfileAlign = File::Spec->catfile($testfileAnalysis, "BOWTIE");
    ok(-e $testfileAlign, "Align directory created");
    my $testfileAlignLog = File::Spec->catfile($testfileAlign, "bowtie_log.txt");
    ok(-e $testfileAlignLog, "Align log file created");
    my $testfileAlignTotal = File::Spec->catfile($testfileAlign, "bowtie.total.mapping.tab");
    ok(-e $testfileAlignTotal, "Align total mapping file created");
    my $testfileAlignQCFile = File::Spec->catfile($testfileQC, "testfile_Bowtie_qc.pdf");
    ok(-e $testfileAlignQCFile, "Align QC pdf file created");
    for my $barcode (split(",", $file->{'Barcodes'})) {
		for my $suffix (qw(unique.output.sort.bam.bai temp.output.sam.gz unique.output.sort.bam unique.output.sam.gz unique.output.bam hitFreq.tab conv.GR.RData chromLoc.tab)) {
    		my $testAlignBarcodeOutputFile = File::Spec->catfile($testfileAlign, "testfile.$barcode.bowtie.$suffix");
    		ok(-e $testAlignBarcodeOutputFile, "Align analysis file testfile.$barcode.bowtie.$suffix created");
    		#ok(-e $testAlignBarcodeOutputFile, "Align analysis file $testAlignBarcodeOutputFile created");
    	}
	}
    
    #features test
    runStep("features", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir --annotationDir=$ANNOTATION");
    my $testfileFeatures = File::Spec->catfile($testfileAnalysis, "miRNA_ANALYSIS");
    ok(-e $testfileFeatures, "Features directory created");
    my $testfileFeatureCount = File::Spec->catfile($testfileFeatures, "testfile.mature.counts.txt");
    ok(-e $testfileFeatureCount, "Align count file created");
}

sub test_pirna_five_prime_barcode_run {
	cleanup();
	prepare();
	
	my $analysisDir = File::Spec->catfile($TEST_RUN, "analysis");
	my @seqfiles = qw(test_file1.txt.gz test_file2.txt.gz test_file3.txt.gz test_file4.txt.gz);
	my $seqLocation = File::Spec->catfile($CWD, "resources", "pirna_5p_barcode");
	for my $seqfile (@seqfiles) {
		my $seqfileLocation = File::Spec->catfile($seqLocation, $seqfile);
		my $newLocation = File::Spec->catfile($TEST_SEQUENCE_STAGING, $seqfile);
		copy($seqfileLocation, $newLocation);
		ok(-e $seqfileLocation, "Test sequence file");
		ok(-e $newLocation, "Test sequence file copy to '$newLocation'");
	}
	
	my $description = [];
	#->{'analysis'} = "pirna_5p_barcode";
	for  (my $i=0; $i <= $#seqfiles; $i++) {
		my $file = {};
		my $seqfile = $seqfiles[$i];
		$file->{"Name"} = "Genome_Research_$i";
		$file->{"File"} = $seqfile;
		$file->{'Geometry'} = "5p_barcode";
		$file->{'Barcodes'} = $i == 0 ? "CTAA" : "ACTA,CTAA";
		$file->{'5p_ad'} = "GTTCAGAGTTCTACAGTCCGACGATC";
		$file->{'3p_ad'} = "ATCTCGTATGCCGTCTTCTGCTTG";
		$file->{'5p_seq_insert'} = "-";
		$file->{'3p_seq_insert'} = "-";
		push(@{$description}, $file);
	}
	my $descriptionFile = writeDescription($description, $TEST_SEQUENCE_STAGING);
	
	my $config = {};
	$config->{'minSize'} = 18;
	$config->{'maxSize'} = 26;
	$config->{'genome'} = "mouse";
	$config->{'ensversion'} = 1;
	$config->{'mismatches'} = 2;
	$config->{'maxHits'} = 20;
	$config->{'sam'} = "FLAG";
	$config->{'feature'} = "miRNA";
	$config->{'mirversion'} = 1;
	$config->{'annot_conflict'} = "merge";
	$config->{'overlap'} = 15;
	my $configFile = writeConfig($config, $TEST_SEQUENCE_STAGING);
    
    #my @analysisNames = map {$_->{"Name"}} @{$description};
    #organise tests
    runStep("organise", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --outDir=$TEST_RUN");
    ok(-e $TEST_ANALYSIS, "Organise directory created");
    for (my $i = 0; $i <= $#{$description}; $i++) {
    	my $name = @{$description}[$i]->{"Name"};
    	my $seqfile = @{$description}[$i]->{"File"};
    	my $testfileAnalysis = File::Spec->catfile($TEST_ANALYSIS, $name);
	    ok(-e $testfileAnalysis, "Organise seqfile subdirectory created");
	    my $testfileMetadata = File::Spec->catfile($testfileAnalysis, "metadata");
	    ok(-e $testfileMetadata, "Organise seqfile metadata subdirectory created");
	    my $testMetadataFile = File::Spec->catfile($testfileMetadata, "metadata.txt");
	    ok(-e $testfileMetadata, "Organise seqfile metadata file created");
	    my $testfileData = File::Spec->catfile($testfileAnalysis, "data");
	    ok(-e $testfileData, "Organise seqfile data subdirectory created");
	    my $testfileDataFile = File::Spec->catfile($testfileData, "analysis_$seqfile");
	    ok(-e $testfileDataFile, "Organise seqfile data file created");	
    }
    
    #reaper tests
    runStep("reaper", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir");
    for (my $i = 0; $i <= $#{$description}; $i++) {
    	my $name = @{$description}[$i]->{"Name"};
    	my $barcodes = @{$description}[$i]->{"Barcodes"};
    	my $seqfile = @{$description}[$i]->{"File"};
    	my $testfileReaper = File::Spec->catfile($TEST_ANALYSIS, $name, "REAPER");
	    ok(-e $testfileReaper, "Reaper analysis directory created");
	    #barcode files
		for my $barcode (split(",", $barcodes)) {
	    	for my $suffix (qw(clean.gz report.input.q report.input.nt report.clean.nt report.clean.len uniquify.log report.clean.trinucl report.clean.annotlen clean.uniquified.fa.gz)) {
	    		my $testReaperBarcodeOutputFile = File::Spec->catfile($testfileReaper, "$name.$barcode.$suffix");
	    		ok(-e $testReaperBarcodeOutputFile, "Reaper file $name.$barcode.$suffix created");
	    	}
	    }
	    #summary files
	    for my $suffix (qw(total.report.miss.nt total.report.miss.len total.report.input.q total.report.input.nt sumstat lint.gz)) {
	    	my $testReaperTotalOutputFile = File::Spec->catfile($testfileReaper, "$name.$suffix");
	    	ok(-e $testReaperTotalOutputFile, "Reaper file $name.$suffix created");
	    }
	    my $testfileQC = File::Spec->catfile($TEST_ANALYSIS, $name, "QC");
	    ok(-e $testfileQC, "Reaper analysis QC directory created");
	    my $testfileReaperQCFile = File::Spec->catfile($testfileQC, $name."_Reaper_qc.pdf");
	    ok(-e $testfileReaperQCFile, "Reaper QC pdf file created");
    }   
    
    #filter tests
    runStep("filter", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir");
    for (my $i = 0; $i <= $#{$description}; $i++) {
    	my $name = @{$description}[$i]->{"Name"};
    	my $barcodes = @{$description}[$i]->{"Barcodes"};
    	my $seqfile = @{$description}[$i]->{"File"};
    	
    	my $testfileFilter = File::Spec->catfile($TEST_ANALYSIS, $name, "PROCESSED");
    	ok(-e $testfileFilter, "Filter directory created");
    
    	my $testfileFilterTrimmed = File::Spec->catfile($testfileFilter, "total.saved.trimmit.tab");
	    ok(-e $testfileFilter, "Filter trimmed file created");
	    my $testfileQC = File::Spec->catfile($TEST_ANALYSIS, $name, "QC");
	    ok(-e $testfileQC, "Reaper analysis QC directory created");
	    my $testfileFilterQC = File::Spec->catfile($testfileQC, $name."_Processed_reads_qc.pdf");
	    ok(-e $testfileFilterQC, "Filter QC file created");
	    
	    #barcode files
		for my $barcode (split(",", $barcodes)) {
			for my $suffix (qw(report.clean.processed.annotlen clean.processed.fa.gz)) {
	    		my $testFilterBarcodeOutputFile = File::Spec->catfile($testfileFilter, "$name.$barcode.$suffix");
	    		ok(-e $testFilterBarcodeOutputFile, "Filter file $name.$barcode.$suffix created");
	    	}
		}
    }    
	
	#align tests
	runStep("align", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir --annotationDir=$ANNOTATION");
    for (my $i = 0; $i <= $#{$description}; $i++) {
    	my $name = @{$description}[$i]->{"Name"};
    	my $barcodes = @{$description}[$i]->{"Barcodes"};
    	my $seqfile = @{$description}[$i]->{"File"};
    	
    	my $testfileAlign = File::Spec->catfile($TEST_ANALYSIS, $name, "BOWTIE");
	    ok(-e $testfileAlign, "Align directory created");
	    my $testfileAlignLog = File::Spec->catfile($testfileAlign, "bowtie_log.txt");
	    ok(-e $testfileAlignLog, "Align log file created");
	    my $testfileAlignTotal = File::Spec->catfile($testfileAlign, "bowtie.total.mapping.tab");
	    ok(-e $testfileAlignTotal, "Align total mapping file created");
	    my $testfileQC = File::Spec->catfile($TEST_ANALYSIS, $name, "QC");
	    ok(-e $testfileQC, "Reaper analysis QC directory created");
	    my $testfileAlignQCFile = File::Spec->catfile($testfileQC, $name."_Bowtie_qc.pdf");
	    ok(-e $testfileAlignQCFile, "Align QC pdf file created");
	    for my $barcode (split(",", $barcodes)) {
			for my $suffix (qw(unique.output.sort.bam.bai temp.output.sam.gz unique.output.sort.bam unique.output.sam.gz unique.output.bam hitFreq.tab conv.GR.RData chromLoc.tab)) {
	    		my $testAlignBarcodeOutputFile = File::Spec->catfile($testfileAlign, "$name.$barcode.bowtie.$suffix");
	    		ok(-e $testAlignBarcodeOutputFile, "Align analysis file $name.$barcode.bowtie.$suffix created");
	    	}
		}
    }
        
    #features test
    runStep("features", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir --annotationDir=$ANNOTATION");
    for (my $i = 0; $i <= $#{$description}; $i++) {
    	my $name = @{$description}[$i]->{"Name"};
    	my $barcodes = @{$description}[$i]->{"Barcodes"};
    	my $seqfile = @{$description}[$i]->{"File"};
    	
    	my $testfileFeatures = File::Spec->catfile($TEST_ANALYSIS, $name, "miRNA_ANALYSIS");
    	ok(-e $testfileFeatures, "Features directory created");
    	my $testfileFeatureCount = File::Spec->catfile($testfileFeatures, "$name.mature.counts.txt");
    	ok(-e $testfileFeatureCount, "Align count file created");
    }
    
}

sub test_small_rna_run {
	cleanup();
	prepare();
	
	my $analysisDir = File::Spec->catfile($TEST_RUN, "analysis");
	my @seqfiles = qw(test_file1.txt.gz test_file3.txt.gz);
	my $seqLocation = File::Spec->catfile($CWD, "resources", "small_rna_no_barcode");
	for my $seqfile (@seqfiles) {
		my $seqfileLocation = File::Spec->catfile($seqLocation, $seqfile);
		my $newLocation = File::Spec->catfile($TEST_SEQUENCE_STAGING, $seqfile);
		copy($seqfileLocation, $newLocation);
		ok(-e $seqfileLocation, "Test sequence file");
		ok(-e $newLocation, "Test sequence file copy to '$newLocation'");
	}
	
	my $description = [];
	#->{'analysis'} = "small_rna_no_barcode";
	for  (my $i=0; $i <= $#seqfiles; $i++) {
		my $file = {};
		my $seqfile = $seqfiles[$i];
		$file->{"Name"} = "Genome_Research_$i";
		$file->{"File"} = $seqfile;
		$file->{'Geometry'} = "no_barcode";
		$file->{'Barcodes'} = "-";
		$file->{'5p_ad'} = "-";
		$file->{'3p_ad'} = "AAAAAAAAAAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA";
		$file->{'5p_seq_insert'} = "-";
		$file->{'3p_seq_insert'} = "-";
		push(@{$description}, $file);
	}
	my $descriptionFile = writeDescription($description, $TEST_SEQUENCE_STAGING);
	
	my $config = {};
	$config->{'minSize'} = 22;
	$config->{'maxSize'} = 32;
	$config->{'genome'} = "human";
	$config->{'ensversion'} = 66;
	$config->{'mismatches'} = 2;
	$config->{'maxHits'} = 20;
	$config->{'sam'} = "FLAG";
	$config->{'feature'} = "miRNA";
	$config->{'mirversion'} = 18;
	$config->{'annot_conflict'} = "merge";
	$config->{'overlap'} = 15;
	my $configFile = writeConfig($config, $TEST_SEQUENCE_STAGING);
    
    #my @analysisNames = map {$_->{"Name"}} @{$description};
    #organise tests
    runStep("organise", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --outDir=$TEST_RUN");
    ok(-e $TEST_ANALYSIS, "Organise directory created");
    for (my $i = 0; $i <= $#{$description}; $i++) {
    	my $name = @{$description}[$i]->{"Name"};
    	my $seqfile = @{$description}[$i]->{"File"};
    	my $testfileAnalysis = File::Spec->catfile($TEST_ANALYSIS, $name);
	    ok(-e $testfileAnalysis, "Organise seqfile subdirectory created");
	    my $testfileMetadata = File::Spec->catfile($testfileAnalysis, "metadata");
	    ok(-e $testfileMetadata, "Organise seqfile metadata subdirectory created");
	    my $testMetadataFile = File::Spec->catfile($testfileMetadata, "metadata.txt");
	    ok(-e $testfileMetadata, "Organise seqfile metadata file created");
	    my $testfileData = File::Spec->catfile($testfileAnalysis, "data");
	    ok(-e $testfileData, "Organise seqfile data subdirectory created");
	    my $testfileDataFile = File::Spec->catfile($testfileData, "analysis_$seqfile");
	    ok(-e $testfileDataFile, "Organise seqfile data file created");	
    }
    
    #reaper tests
    runStep("reaper", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir");
    for (my $i = 0; $i <= $#{$description}; $i++) {
    	my $name = @{$description}[$i]->{"Name"};
    	my $seqfile = @{$description}[$i]->{"File"};
    	my $testfileReaper = File::Spec->catfile($TEST_ANALYSIS, $name, "REAPER");
	    ok(-e $testfileReaper, "Reaper analysis directory created");
	    #barcode files
		for my $suffix (qw(clean.gz report.input.q report.input.nt report.clean.nt report.clean.len uniquify.log report.clean.trinucl report.clean.annotlen clean.uniquified.fa.gz)) {
    		my $testReaperOutputFile = File::Spec->catfile($testfileReaper, "$name.lane.$suffix");
    		ok(-e $testReaperOutputFile, "Reaper file $name.lane.$suffix created");
    	}
    
	    #summary files
	    for my $suffix (qw(sumstat lint.gz)) {
	    	my $testReaperTotalOutputFile = File::Spec->catfile($testfileReaper, "$name.$suffix");
	    	ok(-e $testReaperTotalOutputFile, "Reaper file $name.$suffix created");
	    }
	    my $testfileQC = File::Spec->catfile($TEST_ANALYSIS, $name, "QC");
	    ok(-e $testfileQC, "Reaper analysis QC directory created");
	    my $testfileReaperQCFile = File::Spec->catfile($testfileQC, $name."_Reaper_qc.pdf");
	    ok(-e $testfileReaperQCFile, "Reaper QC pdf file created");
    }   
    
    #filter tests
    runStep("filter", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir");
    for (my $i = 0; $i <= $#{$description}; $i++) {
    	my $name = @{$description}[$i]->{"Name"};
    	my $seqfile = @{$description}[$i]->{"File"};
    	
    	my $testfileFilter = File::Spec->catfile($TEST_ANALYSIS, $name, "PROCESSED");
    	ok(-e $testfileFilter, "Filter directory created");
    
    	my $testfileFilterTrimmed = File::Spec->catfile($testfileFilter, "total.saved.trimmit.tab");
	    ok(-e $testfileFilter, "Filter trimmed file created");
	    my $testfileQC = File::Spec->catfile($TEST_ANALYSIS, $name, "QC");
	    ok(-e $testfileQC, "Reaper analysis QC directory created");
	    my $testfileFilterQC = File::Spec->catfile($testfileQC, $name."_Processed_reads_qc.pdf");
	    ok(-e $testfileFilterQC, "Filter QC file created");
	    
	    #barcode files
		for my $suffix (qw(report.clean.processed.annotlen clean.processed.fa.gz)) {
    		my $testFilterBarcodeOutputFile = File::Spec->catfile($testfileFilter, "$name.lane.$suffix");
    		ok(-e $testFilterBarcodeOutputFile, "Filter file $name.lane.$suffix created");
    	}
    }    
	
	#align tests
	runStep("align", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir --annotationDir=$ANNOTATION");
    for (my $i = 0; $i <= $#{$description}; $i++) {
    	my $name = @{$description}[$i]->{"Name"};
    	my $seqfile = @{$description}[$i]->{"File"};
    	
    	my $testfileAlign = File::Spec->catfile($TEST_ANALYSIS, $name, "BOWTIE");
	    ok(-e $testfileAlign, "Align directory created");
	    my $testfileAlignLog = File::Spec->catfile($testfileAlign, "bowtie_log.txt");
	    ok(-e $testfileAlignLog, "Align log file created");
	    my $testfileAlignTotal = File::Spec->catfile($testfileAlign, "bowtie.total.mapping.tab");
	    ok(-e $testfileAlignTotal, "Align total mapping file created");
	    my $testfileQC = File::Spec->catfile($TEST_ANALYSIS, $name, "QC");
	    ok(-e $testfileQC, "Reaper analysis QC directory created");
	    my $testfileAlignQCFile = File::Spec->catfile($testfileQC, $name."_Bowtie_qc.pdf");
	    ok(-e $testfileAlignQCFile, "Align QC pdf file created");
    	for my $suffix (qw(unique.output.sort.bam.bai temp.output.sam.gz unique.output.sort.bam unique.output.sam.gz unique.output.bam hitFreq.tab conv.GR.RData chromLoc.tab)) {
    		my $testAlignBarcodeOutputFile = File::Spec->catfile($testfileAlign, "$name.lane.bowtie.$suffix");
    		ok(-e $testAlignBarcodeOutputFile, "Align analysis file $name.lane.bowtie.$suffix created");
    	}
    }
        
    #features test
    runStep("features", $configFile, $descriptionFile, "--no-unique --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir --annotationDir=$ANNOTATION");
    for (my $i = 0; $i <= $#{$description}; $i++) {
    	my $name = @{$description}[$i]->{"Name"};
    	my $barcodes = @{$description}[$i]->{"Barcodes"};
    	my $seqfile = @{$description}[$i]->{"File"};
    	
    	my $testfileFeatures = File::Spec->catfile($TEST_ANALYSIS, $name, "miRNA_ANALYSIS");
    	ok(-e $testfileFeatures, "Features directory created");
    	my $testfileFeatureCount = File::Spec->catfile($testfileFeatures, "$name.mature.counts.txt");
    	ok(-e $testfileFeatureCount, "Align count file created");
    }
}

sub test_paired_end_run {
	cleanup();
	prepare();
	
	my $analysisDir = File::Spec->catfile($TEST_RUN, "analysis");
	my @seqfiles = qw(test_file1.fastq.gz test_file2.fastq.gz);
	my $seqLocation = File::Spec->catfile($CWD, "resources", "paired_end", "paired_no_bar");
	for my $seqfile (@seqfiles) {
		my $seqfileLocation = File::Spec->catfile($seqLocation, $seqfile);
		my $newLocation = File::Spec->catfile($TEST_SEQUENCE_STAGING, $seqfile);
		copy($seqfileLocation, $newLocation);
		ok(-e $seqfileLocation, "Test sequence file");
		ok(-e $newLocation, "Test sequence file copy to '$newLocation'");
	}
	
	my $description = [];
	#->{'analysis'} = "small_rna_no_barcode";
	my $file = {};
	$file->{"Name"} = "paired_end";
	$file->{"File"} = join(",", @seqfiles);
	$file->{'Geometry'} = "pair_no_barcode";
	$file->{'Barcodes'} = "-";
	$file->{'5p_ad'} = "TACACTCTTTCCCTACACGACGCTCTTCCGATCT";
	$file->{'3p_ad'} = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
	$file->{'5p_seq_insert'} = "-";
	$file->{'3p_seq_insert'} = "-";
	push(@{$description}, $file);
	my $descriptionFile = writeDescription($description, $TEST_SEQUENCE_STAGING);
	
	my $config = {};
	$config->{'minSize'} = 30;
#	$config->{'maxSize'} = 32;
#	$config->{'five'} = 6;
#	$config->{'genome'} = "human";
#	$config->{'ensversion'} = 66;
#	$config->{'mismatches'} = 2;
#	$config->{'maxHits'} = 20;
	$config->{'fastq'} = "FLAG";
	$config->{'sam'} = "FLAG";
#	$config->{'feature'} = "miRNA";
#	$config->{'mirversion'} = 18;
#	$config->{'annot_conflict'} = "merge";
#	$config->{'overlap'} = 15;
	my $configFile = writeConfig($config, $TEST_SEQUENCE_STAGING);
    
    #my @analysisNames = map {$_->{"Name"}} @{$description};
    #organise tests
    runStep("organise", $configFile, $descriptionFile, "--no-unique --paired --dataDir=$TEST_SEQUENCE_STAGING --outDir=$TEST_RUN");
    ok(-e $TEST_ANALYSIS, "Organise directory created");
    for (my $i = 0; $i <= $#seqfiles; $i++) {
    	my $name = @{$description}[0]->{"Name"}."_".($i+1);
    	my $seqfile = $seqfiles[$i];
    	my $testfileAnalysis = File::Spec->catfile($TEST_ANALYSIS, $name);
	    ok(-e $testfileAnalysis, "Organise seqfile subdirectory created");
	    my $testfileMetadata = File::Spec->catfile($testfileAnalysis, "metadata");
	    ok(-e $testfileMetadata, "Organise seqfile metadata subdirectory created");
	    my $testMetadataFile = File::Spec->catfile($testfileMetadata, "metadata.txt");
	    ok(-e $testfileMetadata, "Organise seqfile metadata file created");
	    my $testfileData = File::Spec->catfile($testfileAnalysis, "data");
	    ok(-e $testfileData, "Organise seqfile data subdirectory created");
	    my $testfileDataFile = File::Spec->catfile($testfileData, "analysis_$seqfile");
	    ok(-e $testfileDataFile, "Organise seqfile data file created");	
    }
    
    #reaper tests
    runStep("reaper", $configFile, $descriptionFile, "--no-unique --paired --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir");
    for (my $i = 0; $i <= $#seqfiles; $i++) {
    	my $name = @{$description}[0]->{"Name"}."_".($i+1);
    	my $seqfile = $seqfiles[$i];
    	my $testfileReaper = File::Spec->catfile($TEST_ANALYSIS, $name, "REAPER");
	    ok(-e $testfileReaper, "Reaper analysis directory created [$testfileReaper]");
	    #barcode files
		for my $suffix (qw(clean.gz report.input.q report.input.nt report.clean.nt report.clean.len report.clean.trinucl)) {
    		my $testReaperOutputFile = File::Spec->catfile($testfileReaper, "$name.lane.$suffix");
    		ok(-e $testReaperOutputFile, "Reaper file $name.lane.$suffix created");
    	}
    
	    #summary files
	    for my $suffix (qw(sumstat lint.gz)) {
	    	my $testReaperTotalOutputFile = File::Spec->catfile($testfileReaper, "$name.$suffix");
	    	ok(-e $testReaperTotalOutputFile, "Reaper file $name.$suffix created");
	    }
	    my $testfileQC = File::Spec->catfile($TEST_ANALYSIS, $name, "QC");
	    ok(-e $testfileQC, "Reaper analysis QC directory created");
	    my $testfileReaperQCFile = File::Spec->catfile($testfileQC, $name."_Reaper_qc.pdf");
	    ok(-e $testfileReaperQCFile, "Reaper QC pdf file created");
    }   
    
    #filter tests
    runStep("filter", $configFile, $descriptionFile, "--no-unique --paired --dataDir=$TEST_SEQUENCE_STAGING --analysisDir=$analysisDir");
    for (my $i = 0; $i <= $#seqfiles; $i++) {
    	my $name = @{$description}[0]->{"Name"}."_".($i+1);
    	my $seqfile = $seqfiles[$i];
     	
    	my $testfileFilter = File::Spec->catfile($TEST_ANALYSIS, $name, "PROCESSED");
    	ok(-e $testfileFilter, "Filter directory created");
    
    	my $testfileFilterTrimmed = File::Spec->catfile($testfileFilter, "total.saved.trimmit.tab");
	    ok(-e $testfileFilter, "Filter trimmed file created");
	    my $testfileQC = File::Spec->catfile($TEST_ANALYSIS, $name, "QC");
	    ok(-e $testfileQC, "Reaper analysis QC directory created");
	    my $testfileFilterQC = File::Spec->catfile($testfileQC, $name."_Processed_reads_qc.pdf");
	    ok(-e $testfileFilterQC, "Filter QC file created");
	    
	    #barcode files
		for my $suffix (qw(tallied.fastq.gz proc.sumstat)) {
    		my $testFilterBarcodeOutputFile = File::Spec->catfile($testfileFilter, "$name.lane.$suffix");
    		ok(-e $testFilterBarcodeOutputFile, "Filter file $name.lane.$suffix created");
    	}
    }    
	
#	#align tests
#	runAlign($TEST_SEQUENCE_STAGING, $TEST_RUN, $configFile, $descriptionFile);
#    for (my $i = 0; $i <= $#{$description}; $i++) {
#    	my $name = @{$description}[$i]->{"Name"};
#    	my $seqfile = @{$description}[$i]->{"File"};
#    	
#    	my $testfileAlign = File::Spec->catfile($TEST_ANALYSIS, $name, "BOWTIE");
#	    ok(-e $testfileAlign, "Align directory created");
#	    my $testfileAlignLog = File::Spec->catfile($testfileAlign, "bowtie_log.txt");
#	    ok(-e $testfileAlignLog, "Align log file created");
#	    my $testfileAlignTotal = File::Spec->catfile($testfileAlign, "bowtie.total.mapping.tab");
#	    ok(-e $testfileAlignTotal, "Align total mapping file created");
#	    my $testfileQC = File::Spec->catfile($TEST_ANALYSIS, $name, "QC");
#	    ok(-e $testfileQC, "Reaper analysis QC directory created");
#	    my $testfileAlignQCFile = File::Spec->catfile($testfileQC, $name."_Bowtie_qc.pdf");
#	    ok(-e $testfileAlignQCFile, "Align QC pdf file created");
#    	for my $suffix (qw(unique.output.sort.bam.bai temp.output.sam.gz unique.output.sort.bam unique.output.sam.gz unique.output.bam hitFreq.tab conv.GR.RData chromLoc.tab)) {
#    		my $testAlignBarcodeOutputFile = File::Spec->catfile($testfileAlign, "$name.lane.bowtie.$suffix");
#    		ok(-e $testAlignBarcodeOutputFile, "Align analysis file $name.lane.bowtie.$suffix created");
#    	}
#    }
#        
#    #features test
#    runFeatures($TEST_SEQUENCE_STAGING, $TEST_RUN, $configFile, $descriptionFile);
#    for (my $i = 0; $i <= $#{$description}; $i++) {
#    	my $name = @{$description}[$i]->{"Name"};
#    	my $barcodes = @{$description}[$i]->{"Barcodes"};
#    	my $seqfile = @{$description}[$i]->{"File"};
#    	
#    	my $testfileFeatures = File::Spec->catfile($TEST_ANALYSIS, $name, "miRNA_ANALYSIS");
#    	ok(-e $testfileFeatures, "Features directory created");
#    	my $testfileFeatureCount = File::Spec->catfile($testfileFeatures, "$name.mature.counts.txt");
#    	ok(-e $testfileFeatureCount, "Align count file created");
#    }
}

#test_Prerequisites();
#test_initialise();
test_five_prime_barcode_run();
#test_pirna_five_prime_barcode_run();
#'test_small_rna_run();
#test_paired_end_run();

#done_testing();
