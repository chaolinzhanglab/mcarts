#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Carp;

use File::Basename;

my $prog = basename ($0);
my $progDir = dirname ($0);

my $verbose = 0;

my $refSpecies = "mm9";
my $trainPosBedFile = "";
my $trainNegBedFile = "";
my $libDataDir = "";

my $existingModel = 0;
my $trainOnly = 0;

my $pattern = "";
my $misMatch = 0;
my $isRegExp = 0;

my $checkMaf = 0;

my $minSiteNum = 3;
my $maxDist = 30;

my $cmdLine = "$prog " . join (" ", @ARGV);

GetOptions ("ref:s"=>\$refSpecies,
			"f:s"=>\$trainPosBedFile,
			"b:s"=>\$trainNegBedFile,
			"lib:s"=>\$libDataDir,
			"w:s"=>\$pattern,
			"m:i"=>\$misMatch,
			"r"=>\$isRegExp,
			"min-site:i"=>\$minSiteNum,
			"max-dist:i"=>\$maxDist,
			"train-only"=>\$trainOnly,
			"exist-model"=>\$existingModel,
			"check-maf"=>\$checkMaf,
			"v"=>\$verbose);

if (@ARGV != 1)
{
	print "search clustered RNA motif sites\n";
	print "Usage: $prog [options] <out dir>\n";
	print "Example1: $prog -v -ref mm9 -f CLIP.pos.bed -b CLIP.neg.bed -lib mm9_mammal_input_data -w YCAY -m 3 --min-site 3 --max-dist 30 out_dir\n";
	print "Example2: $prog -v --exist-model out_dir_from_prev_run\n";

	print "\noptions:\n";
	print " -ref       [string]  : reference species to search ($refSpecies)\n";
	print " -f         [string]  : a BED file specifying positive (foreground) training regions\n";
	print " -b         [string]  : a BED file specifying negative (background) training regions\n";
	print " -lib       [string]  : dir with data library files\n";
	print " -w         [string]  : the consensus motif to search\n";
	print " -m         [int]     : number of mismatches ($misMatch)\n";
	print " -r                   : the motifs provided are regular expressions (will disable -n and -m)\n";
	print " --min-site [int]     : minimum sites in clusters ($minSiteNum)\n";
	print " --max-dist [int]     : max distance allowed in clusters ($maxDist)\n";
	print " --train-only         : train the model only, no prediction\n";
	print " --exist-model        : prediction based on existing model specified in out dir\n";
	print " --check-maf          : check maf files in the library dir\n";
	print " -v                   : verbose\n";
	exit (1);
}

my ($outDir) = @ARGV;

#Carp::croak "failed to create output dir $outDir\n" if -f ($outDir || -d $outDir);

#check integrity of parameters

my ($searchMotif, $getGenicMotif, $formatMotif, $getTrainingSites, $trainHMM) = (1,1,1,1,1);

my $paramFile = "$outDir/params.txt";

if ($existingModel)
{
	Carp::croak "$outDir does not exists\n" unless -d $outDir;
	#Carp::croak "$outDir/BLS does not exists, data inconsistent\n" unless -d "$outDir/BLS";
	Carp::croak "$outDir/formatted does not exists, data inconsistent\n" unless -d "$outDir/formatted";
	Carp::croak "HMM model file $outDir/model.txt does not exist, data inconsistent\n" unless -f "$outDir/model.txt";
	($searchMotif, $getGenicMotif, $formatMotif, $getTrainingSites, $trainHMM) = (0,0,0,0,0);

	#read some command line parameters of the previous run, so that the do not need to be provided again
	my $fin;
	open ($fin, "<$paramFile") || Carp::croak "cannot open file $paramFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\#/;
		next if $line =~/^\s*$/;

		my ($name, $value) = split ("=", $line);
		if ($name eq 'ref')
		{
			$refSpecies = $value;
		}
		elsif ($name eq 'lib')
		{
			$libDataDir = $value;
		}
		elsif ($name eq 'min-site')
		{
			$minSiteNum = $value;
		}
		elsif ($name eq 'max-dist')
		{
			$maxDist = $value;
		}
	}
	close ($fin);
}
else
{
	Carp::croak "No positive training data\n" unless -f $trainPosBedFile;
	Carp::croak "No negative training data\n" unless -f $trainNegBedFile;

	my $ret = system ("mkdir $outDir");
	Carp::croak "failed to create output dir $outDir\n" if $ret != 0;

	$ret = system ("mkdir $outDir/tmp");
	Carp::croak "failed to create output dir $outDir/tmp\n" if $ret != 0;

	my $fout;
	open ($fout, ">$paramFile") || Carp::croak "cannot open file $paramFile to write\n";
	print $fout "#command line: $cmdLine\n";
	print $fout "ref=$refSpecies\n";
	print $fout "lib=$libDataDir\n";
	print $fout "min-site=$minSiteNum\n";
	print $fout "max-dist=$maxDist\n";
	close ($fout);
}

Carp::croak "Library data dir $libDataDir does not exist\n" unless -d $libDataDir;


#----------------------------------------------------------

my $isRegExpFlag = $isRegExp ? "-r" : "";
my $verboseFlag = $verbose ? "-v" : "";

my $motifName = "motif"; 
my $tmpBLSOutDir = "$outDir/BLS";

if ($searchMotif)
{
	print "searching for individual motif sites and evaluate conservation ...\n" if $verbose;
	my $ignoreMafFlag = $checkMaf ? "" : "--ignore-maf";
	my $cmd = "perl $progDir/BLS/searchConservedMotif.pl -name $motifName -ref $refSpecies -w $pattern -m $misMatch $isRegExpFlag $verboseFlag $ignoreMafFlag $libDataDir $tmpBLSOutDir";
	my $ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;
}
else
{
	print "using existing individual motif sites\n" if $verbose;
}

#the output of BLS will be of the form:
#mm9.genic.ext10k.AAAAp.bls.withexonic.txt
#mm9.genic.ext10k.AAAAp.bls.chrom.bed

#split results
#split/bed/mm9.genic.ext10k.maf2fa.AAAAp.mafcoord.0.*.bls.chrom.bed
#split/tab/mm9.genic.ext10k.maf2fa.AAAAp.mafcoord.0.45.bls.withexonic


my $motifWithExonicFile = "$tmpBLSOutDir/$refSpecies.genic.ext10k.$motifName.bls.withexonic.txt";
my $motifBedFile = "$tmpBLSOutDir/$refSpecies.genic.ext10k.$motifName.bls.chrom.bed";

Carp::croak "$motifWithExonicFile does not exists\n" unless -f $motifWithExonicFile;
my $cmd = "wc -l $motifWithExonicFile";
my $numMotifSites = `$cmd`;
chomp $numMotifSites;

if ($numMotifSites =~/^(\d+)\s/)
{
	$numMotifSites = $1;
	print "$numMotifSites individual motif sites found\n";
}
else
{
	Carp::croak "No motif site found\n";
}


#-------------------------------------------------------

my $genicMotifBedFile = "$tmpBLSOutDir/$refSpecies.genic.$motifName.bls.chrom.bed";
if ($getGenicMotif)
{
	print "getting genic motif sites ...\n" if $verbose;
	my $cmd = "perl $progDir/BLS/tagoverlap.pl -big -c $tmpBLSOutDir/cache --keep-score --complete-overlap -region $libDataDir/$refSpecies.genic.bed -ss -d \"#\" $verboseFlag $motifBedFile $genicMotifBedFile.tmp";
	my $ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	$cmd = "sed s/#[a-zA-Z0-9_]*// $genicMotifBedFile.tmp | sort | uniq >$genicMotifBedFile";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	unlink "$genicMotifBedFile.tmp";
}
else
{
	Carp::croak "genic motif site file $genicMotifBedFile does not exist\n" unless -f $genicMotifBedFile;
	print "using existing genic motif sites\n" if $verbose;
}

$cmd = "wc -l $genicMotifBedFile";
$numMotifSites = `$cmd`;
chomp $numMotifSites;

if ($numMotifSites =~/^(\d+)\s/)
{
	$numMotifSites = $1;
	print "$numMotifSites genic motif sites found\n";
}
else
{
	Carp::croak "No motif site found\n";
}


#-----------------------------------------------------------

my $geneContigBedFile = "$libDataDir/$refSpecies.genic.ext10k.bed";
my $geneContigPUDir = "$libDataDir/$refSpecies" . "_genic.ext10k_pu4";
my $tmpFormattedMotifOutDir = "$outDir/formatted";
my $runFormattedMotifOutDir = "$tmpFormattedMotifOutDir/split";
my $formattedMotifFile = "$tmpFormattedMotifOutDir/$refSpecies.genic.ext10k.$motifName.format.txt";

if ($formatMotif)
{
	print "retrieving accessibility information of individual motif sites ...\n" if $verbose;
	my $cmd = "perl $progDir/batchFormatConsMotif.pl $verboseFlag $geneContigBedFile $geneContigPUDir $tmpBLSOutDir/split/tab $tmpFormattedMotifOutDir";
	my $ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	#output file naming: mm9.genic.ext10k.maf2fa.AAAAp.mafcoord.0.*.bls.withexonic.txt
	#

	$cmd = "cat $runFormattedMotifOutDir/$refSpecies.* > $formattedMotifFile";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;
}
else
{
	print "using existing formatted motif sites\n" if $verbose;
}


$cmd = "wc -l $formattedMotifFile";
$numMotifSites = `$cmd`;
chomp $numMotifSites;

if ($numMotifSites =~/^(\d+)\s/)
{
	$numMotifSites = $1;
	print "$numMotifSites formatted motif sites found\n";
}
else
{
	Carp::croak "No motif site found\n";
}





#-----------------------------------------------------------

my $trainPosFile = "$outDir/train_pos.txt";
my $trainNegFile = "$outDir/train_neg.txt";

if ($getTrainingSites)
{
	print "retrieving positive and negative training motif sites ...\n" if $verbose;

	#Pos sites
	my $motif_vs_trainPosBedFile = "$outDir/tmp/motif_vs_train_pos.bed";
	my $cmd = "perl $progDir/BLS/tagoverlap.pl -c $outDir/cache -big -region $trainPosBedFile -ss $verboseFlag -d \"#\" $genicMotifBedFile  $motif_vs_trainPosBedFile";
	my $ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	my $motif_vs_trainPosIdPairFile = "$outDir/tmp/motif_vs_train_pos.idpair";
	$cmd = "awk '{if(\$5==1) {print \$4}}' $motif_vs_trainPosBedFile | awk -F \"#\" '{print \$1\"\\t\"\$2}' | sort | uniq > $motif_vs_trainPosIdPairFile.tmp";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	#remove BLS in the id
	$cmd = "sed 's/\\/\\/[0-9\\.]*\\t/\\t/' $motif_vs_trainPosIdPairFile.tmp > $motif_vs_trainPosIdPairFile";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;
	unlink "$motif_vs_trainPosIdPairFile.tmp";

	my $trainPosFile = "$outDir/train_pos.txt";
	$cmd = "perl $progDir/removeRow.pl -q 3 -r -v $formattedMotifFile $motif_vs_trainPosIdPairFile > $trainPosFile";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	$cmd = "perl $progDir/BLS/selectRow.pl -q 3 $trainPosFile $motif_vs_trainPosIdPairFile > $trainPosFile.tmp";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	$cmd = "paste $trainPosFile.tmp $motif_vs_trainPosIdPairFile | awk '{print \$12\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10}' > $trainPosFile";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;
	unlink "$trainPosFile.tmp";


	####
	#Neg sites
	my $motif_vs_trainNegBedFile = "$outDir/tmp/motif_vs_train_neg.bed";
	$cmd = "perl $progDir/BLS/tagoverlap.pl -c $outDir/cache -big -region $trainNegBedFile -ss $verboseFlag -d \"#\" $genicMotifBedFile  $motif_vs_trainNegBedFile";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	my $motif_vs_trainNegIdPairFile = "$outDir/tmp/motif_vs_train_neg.idpair";
	$cmd = "awk '{if(\$5==1) {print \$4}}' $motif_vs_trainNegBedFile | awk -F \"#\" '{print \$1\"\\t\"\$2}' | sort | uniq > $motif_vs_trainNegIdPairFile.tmp";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	#remove BLS in the id
	$cmd = "sed 's/\\/\\/[0-9\\.]*\\t/\\t/' $motif_vs_trainNegIdPairFile.tmp > $motif_vs_trainNegIdPairFile";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;
	unlink "$motif_vs_trainNegIdPairFile.tmp";

	$cmd = "perl $progDir/removeRow.pl -q 3 -r -v $formattedMotifFile $motif_vs_trainNegIdPairFile > $trainNegFile";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	$cmd = "perl $progDir/BLS/selectRow.pl -q 3 $trainNegFile $motif_vs_trainNegIdPairFile > $trainNegFile.tmp";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	$cmd = "paste $trainNegFile.tmp $motif_vs_trainNegIdPairFile | awk '{print \$12\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10}' > $trainNegFile";
	$ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;
	unlink "$trainNegFile.tmp";
}

#
#
$cmd = "wc -l $trainPosFile";
my $numTrainPosMotifSites = `$cmd`;
chomp $numTrainPosMotifSites;

if ($numTrainPosMotifSites =~/^(\d+)\s/)
{
	$numTrainPosMotifSites = $1;
	print "$numTrainPosMotifSites positive training motif sites found\n";
}
else
{
	Carp::croak "No positive training motif site found\n";
}


$cmd = "wc -l $trainNegFile";
my $numTrainNegMotifSites = `$cmd`;
chomp $numTrainNegMotifSites;

if ($numTrainNegMotifSites =~/^(\d+)\s/)
{
	$numTrainNegMotifSites = $1;
	print "$numTrainNegMotifSites negative training motif sites found\n";
}
else
{
	Carp::croak "No negative training motif site found\n";
}


#------------------------------------------------------------------


my $modelFile = "$outDir/model.txt";
if ($trainHMM)
{
	print "training the HMM ...\n" if $verbose;

	my $cmd = "perl $progDir/HMMMotifCluster.pl $verboseFlag -f $trainPosFile -b $trainNegFile -m $minSiteNum -d $maxDist --save-model $modelFile";
	my $ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;
}
else
{
	Carp::croak "HMM model file $modelFile does not exist\n" unless -f $modelFile;
	print "found existing model file $modelFile\n" if $verbose; 
}

if ($trainOnly)
{
	system ("rm -rf $outDir/tmp");
	exit (0);
}

#------------------------------------------------------------------
print "generating scripts to predict motif clusters ...\n" if $verbose;

my $scriptDir = "$outDir/scripts";
my $ret = 0;

$ret = system ("mkdir $scriptDir") unless -d $scriptDir;
Carp::croak "cannot create dir $scriptDir: $?\n" if $ret != 0;


$cmd = "ls $runFormattedMotifOutDir/$refSpecies.genic.ext10k.maf2fa.$motifName.mafcoord.0.*.bls.withexonic.txt | wc -l";
my $nsplit = `$cmd`;
chomp $nsplit;

my $splitOutDir = "$outDir/out";
$ret = system ("mkdir $splitOutDir") unless -d $splitOutDir;
Carp::croak "cannot create dir $splitOutDir: $?\n" if $ret != 0;

for (my $i = 0; $i < $nsplit; $i++)
{
    print "processing $i ...\n" if $verbose && $i % 100 == 0;

    my $f = "$runFormattedMotifOutDir/$refSpecies.genic.ext10k.maf2fa.$motifName.mafcoord.0.$i.bls.withexonic.txt";
	
    my $outSiteFile = "$splitOutDir/site.$i.bed";
    my $outClusterFile = "$splitOutDir/cluster.$i.bed";

    my $cmd = "perl $progDir/HMMMotifCluster.pl -v --load-model $modelFile -m $minSiteNum $f $outSiteFile $outClusterFile";
    my $scriptFile = "$scriptDir/script.$i.sh";
	my $fout;
	open ($fout, ">$scriptFile");
	print $fout "#!/usr/bin\n";
	print $fout $cmd, "\n";
	close ($fout);
}

#-------------------------------------------------------------------

print "running scripts to predict motif clusters ...\n" if $verbose;

#generate the script list
my $scriptListFile = "$outDir/scripts.list";
system ("ls -rt $scriptDir/*.sh > $scriptListFile");
Carp::croak "The script list file at $scriptListFile was not generated properly\n" unless -f $scriptListFile;

my $qsubCache = "$outDir/qsub";
my $taskName = "mCarts";
$cmd = "rm -rf $qsubCache";
$ret = system ($cmd);
Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

$cmd = "perl $progDir/batchQsub.pl $verboseFlag --wait -c $qsubCache -j $taskName $scriptListFile";
$ret = system ($cmd);
Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;


#------------------------------------------------------------------
print "finishing up ...\n" if $verbose;
my $motifClusterBedFile = "$outDir/cluster.bed";

$cmd = "cat $splitOutDir/cluster.* | grep -v track > $motifClusterBedFile.tmp";
$ret = system ($cmd);
Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

$cmd = "perl $progDir/BLS/contig2genome.pl $verboseFlag $geneContigBedFile $motifClusterBedFile.tmp $motifClusterBedFile";
$ret = system ($cmd);
Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

unlink "$motifClusterBedFile.tmp";
system ("rm -rf $outDir/tmp");



















