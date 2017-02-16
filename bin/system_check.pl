#!/usr/bin/env perl

use Getopt::Long;
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PathUtil qw(findExecutable);

my (@steps) = qw(reaper filter align features complete);
my ($step, $references, $help);

GetOptions(	"step=s" => \$step,
			"references" => \$references,
			"help" => \$help);

if (defined $help) {
	printHelp();
}
#validate arguments
printHelp("No step option supplied. Use 'system_check.pl -s <step>'") unless (defined $step);
my %stepsHash;
@stepsHash{@steps} = ();
printHelp("Unknown step option '$step'") unless (exists $stepsHash{$step});

#check for all steps up to the supplied argument
for(my $i=0; $i < (@steps); $i++) {
	#print "[$i] => $steps[$i]\n";
	if ($i == 0) {
		checkReaper();
		checkR();
		checkRPackages();
	} elsif ($i == 1) {
		checkTally();
	} elsif ($i == 2) {
		checkBowtie();
		checkSamtools();
	}
	if ($steps[$i] eq $step) {
		last;
	}
}

sub checkReaper {
	my $executable = findExecutable("reaper");
	my $exitCode = system("$executable --version");
	die("\nCould not find reaper. To install reaper and for associated documentation please visit http://www.ebi.ac.uk/~stijn/reaper/. You will need to ensure the programme is installed in the file path") if ($exitCode != 0);
	my $version = `$executable --version`;
	print STDERR "\nreaper found. \nVersion: $version\n";
}

sub checkTally {
	my $executable = findExecutable("tally");
	my $exitCode = system("$executable --version");
	die("\nCould not find tally. To install tally and for associated documentation please visit http://www.ebi.ac.uk/~stijn/reaper/.") if ($exitCode != 0);
	my $version = `$executable --version`;
	print STDERR "\ntally found. \nVersion: $version\n";
}

sub checkBowtie {
	my $executable = findExecutable("bowtie");
	my $exitCode = system("$executable --version");
	die("\nCould not find bowtie. Please download and install bowtie (http://bowtie-bio.sourceforge.net/index.shtml).") if ($exitCode != 0);
	my $version = `$executable --version`;
	print STDERR "\nbowtie found. \nVersion: $version\n";
	if ($references) {
		print "Reference:\nUltrafast and memory-efficient alignment of short DNA sequences to the human genome.\nLangmead B., Trapnell C., Pop M. and Salzberg S.\nGenome Biology (2009), 10:R25.\n";
	} 
}

sub checkSamtools {
	my $executable = findExecutable("samtools");
	die("\nCould not find samtools. Please download and install samtools (http://samtools.sourceforge.net/).") 
		unless (-e $executable and -x $executable and -f $executable);
	#my $version = `$executable`;
	print STDERR "\nsamtools found.\n";
	
	if ($references) {
		print "Reference:\nThe Sequence alignment/map (SAM) format and SAMtools.\nLi H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup.\nBioinformatics  (2009), 25:2078-9.\n";
	}
}

sub checkR {
	my $executable = findExecutable("R");
	my $exitCode = system("$executable --version");
	die("\nCould not find R. Please download and install R (http://www.r-project.org/). [$exitCode]") if ($exitCode != 0);
	my $version = `$executable --version`;
	print STDERR "\nR found. \nVersion: $version\n";
}

sub checkRPackages {
	my $executable = findExecutable("R");
	my @packages = qw(GenomicRanges gplots IRanges RColorBrewer ShortRead R.utils);
	foreach my $package (@packages) {
		my $exitCode = system("$executable --vanilla -e 'library(\"$package\")'");
		die("\nCould not find all of the required R libraries. Failed to load $package and any dependencies.") if ($exitCode != 0)
	}
	print STDERR "\nAll the required R libraries are present.\n";
	if ($references) {
		print "References:\n";
        print "R.utils: Various programming utilities. Bengtsson H.\n";
        print "RColorBrewer: ColorBrewer palettes. Neuwirth E.\n";
        print "GenomicRanges: Representation and manipulation of genomic intervals. Aboyoun P., Pages H. and Lawrence M.\n";
        print "gplots: Various R programming tools for plotting data. Warnes G.R.\n";
        print "IRanges: Infrastructure for manipulating intervals on sequences. Pages H., Aboyoun P. and Lawrence M.\n";
        print "ShortRead: a Bioconductor package for input, quality assessment and exploration of high-throughput sequence data.\nMorgan M., Anders S., Lawrence M., Aboyoun P., Pages H. and Gentleman R. Bioinformatics (2009), 25:2607-2608.\n";
	}
}

sub printHelp {
	my ($err) = @_;
	my $exitCode = 0;
	if ($err) {
		print STDERR "$err\n\n" ;
		$exitCode = 1;	
	}
	print "system_check.pl\n";
	print "---------------\n";
	print "Runs system check for a specified seqimp step or all steps if the value is set to'complete'\n";
	print "Arguments:\n";
	print "s: A seqimp step [complete|reaper|filter|align|features\n";
	print "r: Prints references for required executables\n";
	print "h: Prints this message and exit\n";
	print "system_check.pl -s complete:\n";
	exit($exitCode);	
}