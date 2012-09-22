#search motif using a consensus and evaluate branch length score (BLS)
#
# Chaolin Zhang (czhang@rockefeller.edu)
# July 30 2012
#

#!/usr/bin/perl -w
use strict;

use Carp;
use Getopt::Long;
use File::Basename;

use Common;


my $prog = basename ($0);
my $progDir = dirname ($0);

my $pattern = "";
my $isRegExp = 0;
my $numPermute = 0;
my $refSpecies = "mm9";
my $misMatch = 0;
my $searchBothStrand = 0;
my $incSiteLostInRef = 0;
my $motifListFile = "";
my $ignoreMaf = 0; #take whatever maf2fa files that exist, otherwise check maf files and regenerate non-existing maf2fa files

my $taskName = "motif";

my $verbose = 0;

my $cmdLine = "$prog " . join (" ", @ARGV);

GetOptions ("ref:s"=>\$refSpecies,
			"w|word:s"=>\$pattern,
			"r"=>\$isRegExp,
			"m:i"=>\$misMatch,
			"n:i"=>\$numPermute,
			"l:s"=>\$motifListFile,
			"b"=>\$searchBothStrand,
			"inc-loss-in-ref"=>\$incSiteLostInRef,
			"ignore-maf"=>\$ignoreMaf,
			"name:s"=>\$taskName,
			"v"=>\$verbose);


if (@ARGV != 2)
{
	print "search motif using a consensus and evaludate BLS\n";
	print "Usage: $prog [options] <input dir> <out dir>\n";
	print " -ref   [string]  : reference species to search ($refSpecies)\n";
	print " -w     [string]  : the consensus motif to search\n";
	print " -m     [int]     : number of mismatches ($misMatch)\n";
	print " -n     [int]     : the number of permutations of the consensus to search\n";
	print " -l     [string]  : file name with a list of consensus motifs to search (will disable -n)\n";
	print " -r               : the motifs provided are regular expressions (will disable -n and -m)\n"; 
	print " -b               : search both strands\n";
	print " --inc-loss-in-ref: include sites lost in reference genome\n";
	print " --ignore-maf     : ignore maf files and take existing maf2fa files\n";
	print " --name [string]  : task name ($taskName)\n";
	print " -v               : verbose\n";
	exit (1);
}


my ($inputDir, $outputDir) = @ARGV;

Carp::croak "-w and -l cannot be specified at the same time\n" if $pattern && $motifListFile;

$inputDir = getFullPath ($inputDir);


Carp::croak "$outputDir already exists\n" if -d $outputDir;
my $ret = system ("mkdir $outputDir");
Carp::croak "failed to create output dir $outputDir:$?\n" if $ret != 0;

$outputDir = getFullPath ($outputDir);

my $scriptDir = "$outputDir/scripts";
my $ret = system ("mkdir $scriptDir");
Carp::croak "failed to create scripts dir $scriptDir:$?\n" if $ret != 0;



#prepare the list of motifs

print "preparing the list of motif to search ...\n" if $verbose;

my @motifs;
if (-f $motifListFile)
{
	my $fin;
	open ($fin, "<$motifListFile") || Carp::croak "can not open file $motifListFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		push @motifs, $line;
	}
	close ($fin);
	$pattern = $motifs[0];
}
else
{
	push @motifs, $pattern;

	$numPermute = 0 if $isRegExp;

	print "generating $numPermute permuted motifs ...\n" if $numPermute > 0 && $verbose;
	for (my $i = 0; $i < $numPermute; $i++)
	{
		my $c;
		while ($c = `perl $progDir/permuteWord.pl $pattern`)
		{
			chomp $c;
			next if $c eq $pattern;
			last;
		}
		push @motifs, $c;
	}
}

my $nmotifs = @motifs;
print "$nmotifs motifs prepared\n" if $verbose;

my $outMotifFile = "$outputDir/motif.txt";
my $fout;

open ($fout, ">$outMotifFile") || Carp::croak "cannot open file $outMotifFile to write\n";
print $fout "#", $cmdLine, "\n";
print $fout join ("\n", @motifs), "\n";
close ($fout);


print "converting maf files to fasta files (if necessary)...\n" if $verbose;

my @regions = qw(genic.ext10k);
my $nsplit = 0;
my $prefix = $refSpecies;

if ($ignoreMaf)
{
	foreach my $r (@regions)
	{
		my $cmd = "ls $inputDir/$prefix" . "_$r" . "_maf_split/$prefix.$r.*.maf2fa | wc -l";
	
		my $n = `$cmd`;
		chomp $n;
		if ($nsplit == 0)
		{
			$nsplit = $n;
			Carp::croak "cannot detect the number of splits of maf files\n" if $nsplit <= 0;
	
			print "$nsplit splits of maf data detected\n" if $verbose;
		}
		else
		{
			Carp::croak "inconsistency in number of splits for maf files\n" if $n != $nsplit;
		}
	}
}
else
{
	foreach my $r (@regions)
	{
		my $cmd = "ls $inputDir/$prefix" . "_$r" . "_maf_split/$prefix.$r.maf.* | wc -l";
	
		my $n = `$cmd`;
		chomp $n;
		if ($nsplit == 0)
		{
			$nsplit = $n;
			Carp::croak "cannot detect the number of splits of maf files\n" if $nsplit <= 0;
	
			print "$nsplit splits of maf data detected\n" if $verbose;
		}
		else
		{
			Carp::croak "inconsistency in number of splits for maf files\n" if $n != $nsplit;
		}

	
		for (my $i = 0; $i < $nsplit; $i++)
		{
			my $mafFile = "$inputDir/$prefix" . "_$r" . "_maf_split/$prefix.$r.maf.$i";
			my $maf2faFile = "$inputDir/$prefix" . "_$r" . "_maf_split/$prefix.$r.$i.maf2fa";
			next if -f $maf2faFile;
	
			my $cmd = "perl $progDir/maf2fasta.pl  -G -v -s 60 $mafFile $maf2faFile";
			print "$cmd ...\n";
			my $ret = system ($cmd);
			print "CMD $cmd failed: $?\n" if $ret != 0;
		}
	}
}


print "generating scripts for motif search ...\n" if $verbose;

my $speciesFile = "$inputDir/species";
my $treeFile = "$inputDir/tree.nh";

my $exonFile = "$inputDir/$prefix.exon.uniq.bed";
my $fivePrimeUTRFile = "$inputDir/refGene_knownGene.5utr.bed";
my $threePrimeUTRFile = "$inputDir/refGene_knownGene.3utr.bed";


my $resultDir = "$outputDir/split";
my $bedResultDir = "$resultDir/bed";
my $tabResultDir = "$resultDir/tab";

system ("mkdir $resultDir");
system ("mkdir $bedResultDir");
system ("mkdir $tabResultDir");

my $searchBothStrandFlag = $searchBothStrand ? "-b" : ""; 
my $presentInRefSpecies = 1 - $incSiteLostInRef;	#the motif has to be present in reference species

my $exonFile = "$inputDir/$prefix.exon.uniq.bed";
my $fivePrimeUTRFile = "$inputDir/refGene_knownGene.5utr.bed";
my $threePrimeUTRFile = "$inputDir/refGene_knownGene.3utr.bed";

for (my $m = 0; $m < @motifs; $m++)
{
	my $c = $motifs[$m];

	print "$m : $c\n";

	my $runResultDir = "$resultDir/run$m";	
	system ("mkdir $runResultDir");
	
	for (my $i = 0; $i < $nsplit; $i++)
	{
		my $f = "$scriptDir/script$m.$i.sh";
			
		my $fout;
		open ($fout, ">$f") || Carp::croak "can not open file $f to write\n";
		print $fout "#!/bin/sh\n";
		print $fout "\n#split $i\n";

		foreach my $r (@regions)
		{
		
			my $mafFile = "$inputDir/$prefix" . "_$r" . "_maf_split/$prefix.$r.maf.$i";
			my $maf2faFile = "$inputDir/$prefix" . "_$r" . "_maf_split/$prefix.$r.$i.maf2fa";
			
			#search motif
			my $motifFile = "$runResultDir/$prefix.$r.maf2fa.$taskName.$m.$i.bed";
			
			if ($isRegExp)
			{
				print $fout "\nRegExpMatch -c $c $searchBothStrandFlag $maf2faFile > $motifFile\n";
			}
			else
			{
				print $fout "\nPatternMatch -c $c -m $misMatch $searchBothStrandFlag $maf2faFile > $motifFile\n";
			}
			#map genomic coordinates to maf coordinates
			my $motifMafCoordFile = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bed";
			print $fout "perl $progDir/bed2mafcoord.pl -v -trace $motifFile $maf2faFile $motifMafCoordFile\n";
		
			#cluster othologous motifs
			my $motifMafClusterFile = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.tab";
			print $fout "perl $progDir/clusterMafSite.pl -v -species $speciesFile -no-merge-ref -trace -ref $refSpecies $motifMafCoordFile $motifMafClusterFile\n";
		
			#calculate BLS
			my $motifMafClusterBLSFile = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls";
			print $fout "perl $progDir/mafSiteCons.pl -v $treeFile $motifMafClusterFile $motifMafClusterBLSFile\n";

			#filter
			if ($presentInRefSpecies)
			{
				my $motifMafClusterBLSFileTmp = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.tmp";
				print $fout "awk \'{if (\$2>0) {print \$0}}\' $motifMafClusterBLSFile |tail -n +2> $motifMafClusterBLSFileTmp\n";
				print $fout "mv $motifMafClusterBLSFileTmp $motifMafClusterBLSFile\n";
			}

			#bed file with bls
			my $motifMafClusterBLSBedFile = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.bed";
			print $fout "awk \'{print \$1\"//\"\$NF}\' $motifMafClusterBLSFile | awk -F \"//\" \'{print \$1\"\\t\"\$4\"\\t\"\$(NF-1)\"\\t\"\$0\"\\t\"\$NF\"\\t\"\$3}\' > $motifMafClusterBLSBedFile\n";

			#map to genome
			my $seqBedFile = "$inputDir/$prefix.$r.bed";
			my $motifMafClusterBLSGenomeBedFile = "$bedResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.chrom.bed";
			print $fout "perl $progDir/contig2genome.pl $seqBedFile $motifMafClusterBLSBedFile $motifMafClusterBLSGenomeBedFile\n";
			
			#genomic breakdown
			my $motifMafClusterBLSExonBedFile = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.chrom.exon.bed";
			my $motifMafClusterBLS5UTRBedFile = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.chrom.5utr.bed";
			my $motifMafClusterBLS3UTRBedFile = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.chrom.3utr.bed";

			my $ssFlag = $searchBothStrand ? "" : "-ss";
			print $fout "perl $progDir/tagoverlap.pl -c $outputDir/cache $ssFlag -d \"||\" -region $exonFile          -v $motifMafClusterBLSGenomeBedFile $motifMafClusterBLSExonBedFile\n";
			
			my $motifMafClusterBLSExonId = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.chrom.exon.id";
			print $fout "awk \'{print \$4}\' $motifMafClusterBLSExonBedFile | awk -F \'|\' \'{print \$1}\' | awk -F \"//\" '{print \$1\"//\"\$2\"//\"\$3\"//\"\$4\"//\"\$5\"\\t1\"}' | sort | uniq> $motifMafClusterBLSExonId\n";

			print $fout "perl $progDir/tagoverlap.pl -c $outputDir/cache $ssFlag -d \"||\" -region $fivePrimeUTRFile  -v $motifMafClusterBLSGenomeBedFile $motifMafClusterBLS5UTRBedFile\n";
			my $motifMafClusterBLS5UTRId = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.chrom.5utr.id";
			print $fout "awk \'{print \$4}\' $motifMafClusterBLS5UTRBedFile | awk -F \'|\' \'{print \$1}\' | awk -F \"//\" '{print \$1\"//\"\$2\"//\"\$3\"//\"\$4\"//\"\$5\"\\t2\"}' | sort | uniq> $motifMafClusterBLS5UTRId\n";

			print $fout "perl $progDir/tagoverlap.pl -c $outputDir/cache $ssFlag -d \"||\" -region $threePrimeUTRFile -v $motifMafClusterBLSGenomeBedFile $motifMafClusterBLS3UTRBedFile\n";
			my $motifMafClusterBLS3UTRId = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.chrom.3utr.id";
			print $fout "awk \'{print \$4}\' $motifMafClusterBLS3UTRBedFile | awk -F \'|\' \'{print \$1}\' | awk -F \"//\" '{print \$1\"//\"\$2\"//\"\$3\"//\"\$4\"//\"\$5\"\\t3\"}' | sort | uniq> $motifMafClusterBLS3UTRId\n";
			
			#combine the three files
			my $motifMafClusterBLSExonCombineId = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.chrom.exon.combine.id";
			print $fout "cat $motifMafClusterBLS5UTRId $motifMafClusterBLS3UTRId $motifMafClusterBLSExonId >> $motifMafClusterBLSExonCombineId\n";
			
			#remove redundant rows
			print $fout "perl $progDir/uniqRow.pl -c max_num $motifMafClusterBLSExonCombineId $motifMafClusterBLSExonCombineId\n";

			#get intronic or intergenic sites
			my $tmp = "$runResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.chrom.exonic.tmp";
			
			print $fout "perl $progDir/selectRow.pl -p -pt 0 -ss $motifMafClusterBLSExonCombineId $motifMafClusterBLSFile | awk \'{print \$2}\' > $tmp\n";

			my $motifMafClusterBLSExonicFile = "$tabResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.$i.bls.withexonic";
			print $fout "paste $motifMafClusterBLSFile $tmp > $motifMafClusterBLSExonicFile\n";
		}
		print $fout "echo done > $f.done\n";
		close ($fout);
	}
}


print "run jobs ...\n" if $verbose;

my $user = $ENV{'LOGNAME'};

my $sge = 1;
my @qstat = `qstat`;
if ($!)
{
    print "$!\n";
    warn "No SGE service(qsub) in this computer\n" if $sge;
    $sge = 0;
}

#generate the script list
my $scriptListFile = "$outputDir/scripts.list";
system ("ls -rt $scriptDir/*.sh > $scriptListFile");
Carp::croak "The script list file at $scriptListFile was not generated properly\n" unless -f $scriptListFile;

#submit jobs
my $qsubCache = "$outputDir/qsub";
my $queues = "";
system ("perl $progDir/batchQsub.pl -v --wait -c $qsubCache -j $taskName $scriptListFile");


print "combine results ...\n" if $verbose;
for (my $m = 0; $m < @motifs; $m++)
{ 
	foreach my $r (@regions)
	{
		my $splitTag = @motifs > 1 ? "$m." : "";
		my $cmd = "cat $tabResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.*.bls.withexonic > $outputDir/$prefix.$r.$taskName." . $splitTag . "bls.withexonic.txt";
		print "$cmd\n" if $verbose;
		my $ret = system ($cmd);
		Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

		$cmd = "cat $bedResultDir/$prefix.$r.maf2fa.$taskName.mafcoord.$m.*.bls.chrom.bed > $outputDir/$prefix.$r.$taskName." . $splitTag . "bls.chrom.bed";
		print "$cmd\n" if $verbose;
		my $ret = system ($cmd);
		Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;
	}
}



