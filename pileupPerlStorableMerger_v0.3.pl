#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use File::Copy;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use Time::HiRes qw( time );
use Math::Random;
use Math::Round;
use Storable;
use Getopt::Long;
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use Data::Dumper::Names;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a perl script to merge the multiple pileup perl storables generated wiggleToPerlStorable.
#
#	Input
#		--plsIndexListPath=				path; a file contains all paths of the index.hsh.pls, in format of: sample\tplsIndexPath\n;
#		--fastaPath=					reference genome sequence, used for getting the name of the contig, as the reference to look for contig data in the pileup files;
#		--IGVGenomePath=				file path [compulsory]; for coversion of wig into tdf
#		--outDir=						output directory; default = ./wiggleMerger/
#
#		V0.1
#			debut
#
#		v0.2
#			[Sat 14 Sep 2013 19:03:56 CEST] cleaned;
#
#		v0.3
#			[Tue 17 Sep 2013 10:38:17 CEST] Multithread;
#
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-09-17 16:42]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/pileupPerlStorableMerger/v0.3/pileupPerlStorableMerger_v0.3.pl --plsIndexListPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/mergedPileup/mergedSampleList/EHI_Standard.txt --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa --IGVGenomePath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/IGVGenome/EHI_v3.0.genome --outDir=/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/mergedPileup/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/pileupPerlStorableMerger/v0.3/pileupPerlStorableMerger_v0.3.pl
#	--plsIndexListPath=/Volumes/C_Analysis/NGS/results/E014_heatShock_polyA_basic_min35nt/T8H_mergeList.txt
#	--fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa
#	--IGVGenomePath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/IGVGenome/EHI_v3.0.genome
#	--outDir=/Volumes/C_Analysis/NGS/results/E014_heatShock_polyA_basic_min35nt/mergedPileup/T8H/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $ARGVStr = join "\n", (&currentTime(), abs_path($0), @ARGV);#->245
my $globalScriptDirPath = dirname(rel2abs($0));
open DEBUGLOG, ">", "$globalScriptDirPath/debug.log.txt";
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|414, readParameters|595
#	secondaryDependOnSub: currentTime|245
#
#<section ID="startingTasks" num="0">
&printCMDLogOrFinishMessage("CMDLog");#->414
#----------Read parameters ----------#
my ($plsIndexListPath, $fastaPath, $IGVGenomePath, $outDir) = &readParameters();#->595
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="1">
my ($plsIndexListName, undef, undef) = fileparse($plsIndexListPath, qr/\.[^.]*/);
my $maxThread = 12;
my $paramTag = "$plsIndexListName";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="2">
my @mkDirAry;
my $resultDir = "$outDir/$paramTag/"; push @mkDirAry, $resultDir;
my $resultStorableDir = "$resultDir/storable/"; push @mkDirAry, $resultStorableDir;
my $resultLogDir = "$resultDir/log/"; push @mkDirAry, $resultLogDir;
my $resultWigDir = "$resultDir/wig/"; push @mkDirAry, $resultWigDir;
foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineOutFilePath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutFilePath" num="3">
my $calledCMDPath = "$resultDir/called.CMD.txt";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_processInputData
#	primaryDependOnSub: getCntgCovPlsPaths|290, readContigLengthFromFasta|538
#	secondaryDependOnSub: reportStatus|628
#
#<section ID="processInputData" num="4">
#----------Read contig names
my ($cntgLenHsh_ref) = &readContigLengthFromFasta($fastaPath);#->538
#----------Check pls index
my ($cntgPosDataPlsPathBySamByCntgHsh_ref) = &getCntgCovPlsPaths($plsIndexListPath, $outDir);#->290
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_mergeContigs
#	primaryDependOnSub: mergeCntgCovPls|330
#	secondaryDependOnSub: checkRunningThreadAndWaitToJoin|202, generateThreadHshWithRandomItem|263, reportStatus|628
#
#<section ID="mergeContigs" num="5">
my ($mergedCntgPosDataPlsPathHsh_ref) = &mergeCntgCovPls($cntgPosDataPlsPathBySamByCntgHsh_ref, $cntgLenHsh_ref, $resultStorableDir, $maxThread);#->330
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_outputWig
#	primaryDependOnSub: printWigAndTDFFromCntgPls|466
#	secondaryDependOnSub: checkRunningThreadAndWaitToJoin|202, convertWigToTDF|226, printWiggleSingleTrackFromCntgCovPlsPathHsh|496, reportStatus|628
#
#<section ID="outputWig" num="6">
#----print merged wigs
&printWigAndTDFFromCntgPls($mergedCntgPosDataPlsPathHsh_ref, $IGVGenomePath, $resultWigDir);#->466
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 7_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|414, printCalledCMD|447
#	secondaryDependOnSub: currentTime|245
#
#<section ID="finishingTasks" num="7">
&printCalledCMD($ARGVStr, $calledCMDPath);#->447
&printCMDLogOrFinishMessage("finishMessage");#->414
close DEBUGLOG;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	fasta [n=1]:
#		readContigLengthFromFasta
#
#	general [n=5]:
#		currentTime, printCMDLogOrFinishMessage, printCalledCMD
#		readParameters, reportStatus
#
#	multithread [n=2]:
#		checkRunningThreadAndWaitToJoin, generateThreadHshWithRandomItem
#
#	reporting [n=1]:
#		currentTime
#
#	specific [n=1]:
#		mergeCntgCovPls
#
#	storable [n=2]:
#		getCntgCovPlsPaths, mergeCntgCovPls
#
#	wiggle [n=3]:
#		convertWigToTDF, printWigAndTDFFromCntgPls, printWiggleSingleTrackFromCntgCovPlsPathHsh
#
#====================================================================================================================================================#

sub checkRunningThreadAndWaitToJoin {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: >none
#	appearInSub: mergeCntgCovPls|330, printWigAndTDFFromCntgPls|466
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_mergeContigs|138, 6_outputWig|148
#	input: none
#	output: none
#	toCall: &checkRunningThreadAndWaitToJoin();
#	calledInLine: 409, 487, 492
#....................................................................................................................................................#

	my @runningThrAry = threads->list(threads::running);
	my @joinableThrAry = threads->list(threads::joinable);
	while (@runningThrAry or @joinableThrAry) {
		@runningThrAry = threads->list(threads::running);
		@joinableThrAry = threads->list(threads::joinable);
		foreach my $joinableThr (@joinableThrAry) {
			$joinableThr->detach() if not $joinableThr->is_running();
		}
		sleep 2;
	}
}
sub convertWigToTDF {
#....................................................................................................................................................#
#	subroutineCategory: wiggle
#	dependOnSub: >none
#	appearInSub: printWigAndTDFFromCntgPls|466
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_outputWig|148
#	input: $IGVGenomePath, $wigPath
#	output: 
#	toCall: &convertWigToTDF($wigPath, $IGVGenomePath);
#	calledInLine: 490, 491
#....................................................................................................................................................#
	my ($wigPath, $IGVGenomePath) = @_;

	my ($wigName, $wigDir, undef) = fileparse($wigPath, qr/\.[^.]*/);
	system ("igvtools toTDF $wigPath $wigDir$wigName.tdf $IGVGenomePath 2&>/dev/null;");

	return ();
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: >none
#	appearInSub: printCMDLogOrFinishMessage|414, reportStatus|628
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|76, 7_finishingTasks|159
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 65, 434, 437, 442, 644
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub generateThreadHshWithRandomItem {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: >none
#	appearInSub: mergeCntgCovPls|330
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_mergeContigs|138
#	input: $itemAry_ref, $maxThread
#	output: $randItemForThrHsh_ref
#	toCall: my ($randItemForThrHsh_ref) = &generateThreadHshWithRandomItem($maxThread, $itemAry_ref);
#	calledInLine: 363
#....................................................................................................................................................#

	my ($maxThread, $itemAry_ref) = @_;

	my @shuffleItemAry = shuffle(@{$itemAry_ref});
	my $threadNum = 1;
	my $randItemForThrHsh_ref = {};
	foreach my $item (@shuffleItemAry) {
		$threadNum = 1 if $threadNum > $maxThread;
		push @{$randItemForThrHsh_ref->{$threadNum}}, $item;
		$threadNum++;
	}
	
	return ($randItemForThrHsh_ref);

}
sub getCntgCovPlsPaths {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: reportStatus|628
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|125
#	secondaryAppearInSection: >none
#	input: $outDir, $plsIndexListPath
#	output: $cntgPosDataPlsPathBySamByCntgHsh_ref
#	toCall: my ($cntgPosDataPlsPathBySamByCntgHsh_ref) = &getCntgCovPlsPaths($plsIndexListPath, $outDir);
#	calledInLine: 133
#....................................................................................................................................................#
	
	my ($plsIndexListPath, $outDir) = @_;
	
	my %plsIndexPathHsh = ();
	my @sampleNameAry = ();
	my $cntgPosDataPlsPathBySamByCntgHsh_ref = {};
	
	open (IDXLISTPATH, $plsIndexListPath);
	while (<IDXLISTPATH>) {
		chomp;
		next if $_ =~ m/^#/;
		my ($sample, $plsIndexPath) = split /\t/;
		my ($plsIndexName, $plsIndexDir, $plsIndexSuffix) = fileparse($plsIndexPath, qr/\.[^.]*/);
		system ("gzip -d -f $plsIndexPath.gz") if (-s "$plsIndexPath.gz");
		die "wiggle file of $sample doesn't exists\n" unless (-s $plsIndexPath);
		&reportStatus("$sample .pls index checked", 10, "\n");#->628
		$plsIndexPathHsh{$sample} = $plsIndexPath;
		push @sampleNameAry, $sample;
		my %plsIndexHsh = %{retrieve($plsIndexPath)};
		my (undef, $cntgCovStroableDir, undef) = fileparse($plsIndexPath, qr/\.[^.]*/);
		foreach my $cntg (keys %plsIndexHsh) {
			$cntgPosDataPlsPathBySamByCntgHsh_ref->{$sample}{$cntg} = "$cntgCovStroableDir/$plsIndexHsh{$cntg}";
		}
	}
	close IDXLISTPATH;

	return ($cntgPosDataPlsPathBySamByCntgHsh_ref);
}
sub mergeCntgCovPls {
#....................................................................................................................................................#
#	subroutineCategory: storable, specific
#	dependOnSub: checkRunningThreadAndWaitToJoin|202, generateThreadHshWithRandomItem|263, reportStatus|628
#	appearInSub: >none
#	primaryAppearInSection: 5_mergeContigs|138
#	secondaryAppearInSection: >none
#	input: $cntgLenHsh_ref, $cntgPosDataPlsPathBySamByCntgHsh_ref, $maxThread, $resultStorableDir
#	output: $mergedCntgPosDataPlsPathHsh_ref
#	toCall: my ($mergedCntgPosDataPlsPathHsh_ref) = &mergeCntgCovPls($cntgPosDataPlsPathBySamByCntgHsh_ref, $cntgLenHsh_ref, $resultStorableDir, $maxThread);
#	calledInLine: 143
#....................................................................................................................................................#

	my ($cntgPosDataPlsPathBySamByCntgHsh_ref, $cntgLenHsh_ref, $resultStorableDir, $maxThread) = @_;
	
	my $mergedCntgPosDataPlsPathHsh_ref = ();

	#----create empty ary
	my $cntgCovStroableDir = "$resultStorableDir/merged";
	system "mkdir -pm 777 $cntgCovStroableDir";
	my $mergedPlsIndexHsh_ref = {};
	
	foreach my $cntg (keys %{$cntgLenHsh_ref}) {
		&reportStatus("Creating empty storabe for $cntg", 10, "\r");#->628
		my $cntgCovAry_ref = [];
		push @{$cntgCovAry_ref}, undef foreach (1..$cntgLenHsh_ref->{$cntg});
		my $cntgCovStroableName = "$cntg.ary.pls";
		my $cntgCovStroablePath = "$cntgCovStroableDir/$cntgCovStroableName";
		$mergedCntgPosDataPlsPathHsh_ref->{$cntg} = $cntgCovStroablePath;
		$mergedPlsIndexHsh_ref->{$cntg} = $cntgCovStroableName;
		store($cntgCovAry_ref, "$cntgCovStroablePath");
	}
	store($mergedPlsIndexHsh_ref, "$cntgCovStroableDir/index.hsh.pls");

	my ($randCntgInThrHsh_ref) = &generateThreadHshWithRandomItem($maxThread, [keys %{$cntgLenHsh_ref}]);#->263
	my $cntgProc :shared = 0;
	my %threadHsh = ();

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThrHsh_ref}) {
		my $cntgAry_ref = $randCntgInThrHsh_ref->{$threadNum};
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
			sub {

				my ($cntgAry_ref) = @_;

				foreach my $cntg (@{$cntgAry_ref}) {
					$cntgProc++;
					&reportStatus("$cntgProc processed", 30, "\r");#->628
					my $mergedCntgCovPlsPath = $mergedCntgPosDataPlsPathHsh_ref->{$cntg};
					system ("gzip -d -f $mergedCntgCovPlsPath.gz") if (-s "$mergedCntgCovPlsPath.gz");
					my $mergeCntgCovAry_ref = retrieve($mergedCntgCovPlsPath);
					foreach my $sample (keys %{$cntgPosDataPlsPathBySamByCntgHsh_ref}) {
						my $cntgPosDataPlsPath = $cntgPosDataPlsPathBySamByCntgHsh_ref->{$sample}{$cntg};
						system ("gzip -d -f $cntgPosDataPlsPath.gz") if (-s "$cntgPosDataPlsPath.gz");
						my $cntgCovAry_ref = retrieve($cntgPosDataPlsPath);

						my (undef, $cntgCovStroableDir, undef) = fileparse($cntgPosDataPlsPath, qr/\.[^.]*/);

						foreach my $index (0..$#{$mergeCntgCovAry_ref}) {
							next if (not defined $cntgCovAry_ref->[$index]);
							my ($mergedPlusCov, $mergedMinusCov); 
							if (defined $mergeCntgCovAry_ref->[$index]) {
								my($basePlusCov, $baseMinusCov) = split /,/, $mergeCntgCovAry_ref->[$index];
								my ($samplePlusCov, $sampleMinusCov) = split /,/, $cntgCovAry_ref->[$index];
								$mergedPlusCov = $samplePlusCov + $basePlusCov;
								$mergedMinusCov = $sampleMinusCov + $baseMinusCov;
							} else {
								($mergedPlusCov, $mergedMinusCov) = split /,/, $cntgCovAry_ref->[$index];
							}

							$mergeCntgCovAry_ref->[$index] = join ",", ($mergedPlusCov, $mergedMinusCov);
						}
					}
					store($mergeCntgCovAry_ref, $mergedCntgCovPlsPath);
				}
			}#---end of sub
			,($cntgAry_ref)
		); #---end of threads->new
	}#---end of foreach my $threadNum

	&checkRunningThreadAndWaitToJoin();#->202
	
	return ($mergedCntgPosDataPlsPathHsh_ref);
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|245
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|76, 7_finishingTasks|159
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 81, 165
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->245
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->245
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->245
		print "=========================================================================\n\n";
	}
}
sub printCalledCMD {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 7_finishingTasks|159
#	secondaryAppearInSection: >none
#	input: $ARGVStr, $calledCMDPath
#	output: 
#	toCall: &printCalledCMD($ARGVStr, $calledCMDPath);
#	calledInLine: 164
#....................................................................................................................................................#
	my ($ARGVStr, $calledCMDPath) = @_;
	
	open CALLEDCMD, ">", $calledCMDPath;
	print CALLEDCMD join "", ($ARGVStr), "\n";
	close CALLEDCMD;
	return ();
}
sub printWigAndTDFFromCntgPls {
#....................................................................................................................................................#
#	subroutineCategory: wiggle
#	dependOnSub: checkRunningThreadAndWaitToJoin|202, convertWigToTDF|226, printWiggleSingleTrackFromCntgCovPlsPathHsh|496, reportStatus|628
#	appearInSub: >none
#	primaryAppearInSection: 6_outputWig|148
#	secondaryAppearInSection: >none
#	input: $IGVGenomePath, $cntgPosDataPlsPathHsh_ref, $resultWigDir
#	output: none
#	toCall: &printWigAndTDFFromCntgPls($cntgPosDataPlsPathHsh_ref, $IGVGenomePath, $resultWigDir);
#	calledInLine: 154
#....................................................................................................................................................#

	my ($cntgPosDataPlsPathHsh_ref, $IGVGenomePath, $resultWigDir) = @_;

	my $zipPls = "no";
	my $zipWig = "yes";
	my $plusWigPath = "$resultWigDir/merged.plus.wig.gz";
	my $minusWigPath = "$resultWigDir/merged.minus.wig.gz";
	&reportStatus("Issued 2 threads to print plus and minue merged wiggle", 10, "\n");#->628
	threads->new(\&printWiggleSingleTrackFromCntgCovPlsPathHsh, ($cntgPosDataPlsPathHsh_ref, 0, $plusWigPath, $zipPls, $zipWig));#->496
	threads->new(\&printWiggleSingleTrackFromCntgCovPlsPathHsh, ($cntgPosDataPlsPathHsh_ref, 1, $minusWigPath, $zipPls, $zipWig));#->496
	&checkRunningThreadAndWaitToJoin();#->202
	
	&reportStatus("Issued 2 threads to print convert wiggle to tdf", 10, "\n");#->628
	threads->new(\&convertWigToTDF, ($plusWigPath, $IGVGenomePath));#->226
	threads->new(\&convertWigToTDF, ($minusWigPath, $IGVGenomePath));#->226
	&checkRunningThreadAndWaitToJoin();#->202

}
sub printWiggleSingleTrackFromCntgCovPlsPathHsh {
#....................................................................................................................................................#
#	subroutineCategory: wiggle
#	dependOnSub: >none
#	appearInSub: printWigAndTDFFromCntgPls|466
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_outputWig|148
#	input: $aryIndex, $cntgPosDataPlsPathHsh_ref, $wigPath, $zipPls, $zipWig
#	output: none
#	toCall: &printWiggleSingleTrackFromCntgCovPlsPathHsh($cntgPosDataPlsPathHsh_ref, $aryIndex, $wigPath, $zipPls, $zipWig);
#	calledInLine: 485, 486
#....................................................................................................................................................#
	
	my ($cntgPosDataPlsPathHsh_ref, $aryIndex, $wigPath, $zipPls, $zipWig) = @_;
	
	if ($zipWig eq 'no') {
		open (WIGGLE, ">", $wigPath);
	} elsif ($zipWig eq 'yes') {
		open (WIGGLE, "| gzip -fc >$wigPath");
	}
	
	foreach my $cntg (sort keys %{$cntgPosDataPlsPathHsh_ref}) {

		print WIGGLE "variableStep chrom=$cntg span=1\n";

		my $cntgPosDataPlsPath = "$cntgPosDataPlsPathHsh_ref->{$cntg}";
 		system ("gzip -df $cntgPosDataPlsPath.gz") if (-s "$cntgPosDataPlsPath.gz" and $zipPls eq 'yes');
		my $cntgCovAry_ref = retrieve($cntgPosDataPlsPath);
		system ("gzip -f $cntgPosDataPlsPath") if (-s $cntgPosDataPlsPath and $zipPls eq 'yes');
		for my $i (0..$#{$cntgCovAry_ref}) {
			if ($cntgCovAry_ref->[$i]) {
				my @tmpCovAry = split /,/, $cntgCovAry_ref->[$i];
				my $cov = $tmpCovAry[$aryIndex];
				if ($cov > 0) {
					my $pos = $i + 1;
					print WIGGLE join '', ((join "\t", ($pos, $cov)), "\n");
				}
			}
		}
	}
	close WIGGLE;
}
sub readContigLengthFromFasta {
#....................................................................................................................................................#
#	subroutineCategory: fasta
#	dependOnSub: reportStatus|628
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|125
#	secondaryAppearInSection: >none
#	input: $fastaPath
#	output: \%cntgLenHsh
#	toCall: my (\%cntgLenHsh) = &readContigLengthFromFasta($fastaPath);
#	calledInLine: 131
#....................................................................................................................................................#

	my ($fastaPath) = @_;
	
	&reportStatus("Reading $fastaPath for contig names", 10, "\n");#->628
	open (INFILE, $fastaPath);
	my (%cntgLenHsh, $seqName, $length, $seq);

	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
			
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			$seqName = $curntLine;
			$seqName =~ s/>//g; #---remove space
			$cntgLenHsh{$seqName} = 0;

		} else {#--seq line
			$cntgLenHsh{$seqName} = $cntgLenHsh{$seqName} + length ($curntLine);
		}
			
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seqName = $nextLine;
			$seqName =~ s/>//g; #---remove space
			$cntgLenHsh{$seqName} = 0;

		} elsif (eof(INFILE)) {#---this is the last line
			$cntgLenHsh{$seqName} = $cntgLenHsh{$seqName} + length ($nextLine);
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}
	close INFILE;

	my $contigNum = keys %cntgLenHsh;

	&reportStatus("Totally $contigNum contig names stored", 10, "\n");#->628
	
	#---print to check
	#foreach my $cntg (sort {$a cmp $b} keys %cntgLenHsh) {print $cntg."\t".$cntgLenHsh{$cntg}."\n";}

	return (\%cntgLenHsh);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|76
#	secondaryAppearInSection: >none
#	input: none
#	output: $IGVGenomePath, $fastaPath, $outDir, $plsIndexListPath
#	toCall: my ($plsIndexListPath, $fastaPath, $IGVGenomePath, $outDir) = &readParameters();
#	calledInLine: 83
#....................................................................................................................................................#
	my ($plsIndexListPath, $fastaPath, $IGVGenomePath, $outDir);
	
	my $dirPath = dirname(rel2abs($0));
	$outDir = "$dirPath/pileupPerlStorableMerger/";

	GetOptions 	("plsIndexListPath=s"  => \$plsIndexListPath,
				 "fastaPath=s"  => \$fastaPath,
				 "IGVGenomePath=s"  => \$IGVGenomePath,
				 "outDir:s"  => \$outDir)

	or die		("Error in command line arguments\n");

	#---check the files
	foreach my $fileToCheck ($plsIndexListPath, $fastaPath, $IGVGenomePath) {
		die "Can't read $fileToCheck" if not -s $fileToCheck;
	}
	
	chop $outDir if ($outDir =~ m/\/$/); #---remove the last slash
	
	return ($plsIndexListPath, $fastaPath, $IGVGenomePath, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|245
#	appearInSub: getCntgCovPlsPaths|290, mergeCntgCovPls|330, printWigAndTDFFromCntgPls|466, readContigLengthFromFasta|538
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_processInputData|125, 5_mergeContigs|138, 6_outputWig|148
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 316, 352, 376, 484, 489, 552, 587
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->245

	return ();
}

exit;
