#!/usr/bin/perl -w

use strict;
use Carp;
use File::Basename;
use Getopt::Long;

my $prog = basename ($0);
my $progDir = dirname ($0);
my $ext = 0;
my $noCoord = 0;
my $endToTruncate = 0;

my $verbose = 0;

GetOptions ("e:i"=>\$ext,
	"k"=>\$noCoord,
	"t:i"=>\$endToTruncate,
	"v"=>\$verbose);

if (@ARGV != 4)
{
	print "format motif sites for mCarts\n";
	print "$prog [options] <seq.ref.bed> <seq_ref_unpaired_dir> <site_bls_withexonic_dir> <out_dir>\n";
	print " -ext    [int]: extension on each site of the exon (used to define uif, exon, dif) ($ext)\n";
	print " -i      [int]: number of nucleotide to truncate seq in PU files ($endToTruncate)\n";
	print " --no-coord   : do not keep coordinate information\n";
	print " -v           : verbose\n";
	exit (1);
}

my ($seqRefGenomeBedFile, $seqRefProbUnpairDir, $siteBLSWithExonicDir, $outDir) = @ARGV;

Carp::croak "dir does not exists\n" unless (-d $seqRefProbUnpairDir) && (-d $siteBLSWithExonicDir) && (-f $seqRefGenomeBedFile); 

my $ret = system ("mkdir $outDir");
Carp::croak "cannot create dir $outDir : $?\n" if $ret != 0;

my @siteFiles = `ls $siteBLSWithExonicDir`;
my @probUnpairFiles = `ls $seqRefProbUnpairDir`;

chop @siteFiles;
chop @probUnpairFiles;

Carp::croak "the number of split in $seqRefProbUnpairDir and $siteBLSWithExonicDir are different\n" unless $#siteFiles == $#probUnpairFiles;

my $nsplit = @siteFiles;

print "$nsplit splits detected\n" if $verbose;

my $f = $siteFiles[0];

$f=~/^(.*?mafcoord\.\d+)\.\d+\.(.*?)$/;

my $siteFilePrefix = $1;
my $siteFileSuffix = $2;

$f = $probUnpairFiles[0];
$f=~/^(.*?\.pu)\.\d+$/;

my $probUnpairFilePrefix = $1;

print "site prefix=$siteFilePrefix, suffix=$siteFileSuffix, pu prefix=$probUnpairFilePrefix\n" if $verbose;

print "generating scripts ...\n" if $verbose;

my $scriptDir = "$outDir/scripts";
$ret = system ("mkdir $scriptDir");
print "cannot create dir $scriptDir: $?\n" if $ret != 0;

my $runOutDir = "$outDir/split";
$ret = system ("mkdir $runOutDir");
print "cannot create dir $runOutDir: $?\n" if $ret != 0;


my $i = 0;

my $verboseFlag = $verbose ? "-v" : "";

foreach my $siteFile (@siteFiles)
{
	print "processing $i ...\n" if $verbose && $i % 100 == 0;
	$i++;

	$siteFile =~/^.*?mafcoord\.\d+\.(\d+)\..*?$/;
	my $iter = $1;

	$siteFile = "$siteBLSWithExonicDir/$siteFile";

	my $probUnpairFile = "$seqRefProbUnpairDir/" . $probUnpairFilePrefix . "." . $iter;
	my $outFile = "$runOutDir/$siteFilePrefix.$iter.$siteFileSuffix.txt";
	my $outUIFFile = "$runOutDir/$siteFilePrefix.$iter.$siteFileSuffix.uif1k.txt";
	my $outExonFile = "$runOutDir/$siteFilePrefix.$iter.$siteFileSuffix.exon.txt";
	my $outDIFFile = "$runOutDir/$siteFilePrefix.$iter.$siteFileSuffix.dif1k.txt";

	my $keepCoordFlag = $noCoord ? '' : "-keepcoord";
	my $cmd = "perl ~/scripts/formatConsMotif.pl -ext $ext $keepCoordFlag -truncate $endToTruncate $verboseFlag $seqRefGenomeBedFile $probUnpairFile $siteFile $outFile";

	my $scriptFile = "$scriptDir/script.$iter.sh";
	my $fout;
    open ($fout, ">$scriptFile");
    print $fout "#!/usr/bin\n";
    print $fout $cmd, "\n";

	#these are obsolete code, but we keep it here for now
	if ($noCoord)
	{
		print $fout "grep uif $outFile > $outUIFFile\n";
		print $fout "grep exon $outFile > $outExonFile\n";
		print $fout "grep dif $outFile > $outDIFFile\n";
		print $fout "rm $outFile\n";
	}
	close ($fout);
}

#generate the script list
my $scriptListFile = "$outDir/scripts.list";
system ("ls -rt $scriptDir/*.sh > $scriptListFile");
Carp::croak "The script list file at $scriptListFile was not generated properly\n" unless -f $scriptListFile;

my $qsubCache = "$outDir/qsub";
my $taskName = "FormatMotif";

my $cmd = "perl $progDir/batchQsub.pl $verboseFlag --wait -c $qsubCache -j $taskName $scriptListFile";
$ret = system ($cmd);
Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

#system ("rm -rf $scriptDir");
#system ("rm -rf $qsubCache");
#unlink $scriptListFile



