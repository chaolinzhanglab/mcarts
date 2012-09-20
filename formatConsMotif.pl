#!/usr/bin/perl -w

#use the output of the BLS pipeline to generate file format required by HMMConsMotif.pl

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Carp;

use Common;
use Bed;

my $verbose = 0;
my $ext = 0;
my $prog = basename ($0);
my $keepCoord = 0;	#if 1, $ext is not required
my $puWin = 0;
my $endToTruncate = 0;

my $cmd = join (" ", $prog, @ARGV);


GetOptions ("ext:i"=>\$ext,
	"keepcoord"=>\$keepCoord,
	"win:i"=>\$puWin,
	"truncate:i"=>\$endToTruncate,
	"v|verbose"=>\$verbose);

if (@ARGV != 4)
{
	print "generate formated motif file\n";
	print "Usage: $prog [options] <seq.ref.bed (in)> <seq.ref.unpaired.txt> <site_bls_withexonic.txt (in)> <site.txt (out)>\n";
	print "options:\n";
	print " -ext   [int]   : extension on each site of the exon (used to define uif, exon, dif) ($ext)\n";
	print " -keepcoord     : keep the original coordinates\n";
	print " -win   [int]   : window size to estimate pu (obsolete)\n";
	print " -truncate [int]: end to truncate (n=$endToTruncate)\n";
	print " -v             : verbose\n";
	exit (0);
}

Carp::croak "the option -win is obosolete, use -truncate instead\n" if ($puWin != 0);

print "CMD= $cmd\n" if $verbose;

my %seqToRemove = (); #"chr13_f_c347"=> 1, "chr3_f_c395"=>1, "chr6_f_c467"=>1, "chrM_f_c0"=>1, "chr13_r_c364"=>1, "chr13_random_r_c2"=> 1, "chr5_random_r_c0"=>1, "chrUn_random_r_c0"=>1, "chrY_r_c0"=>1); #, "uc009vev.1_0"=>1, NM_001039244_2=>1, "uc007rzm.1_1"=>1); 

#hacker: these sequences should be removed because they are close to the termini of chromosomes and the extension is incorrect


my ($seqRefGenomeBedFile, $seqRefProbUnpairFile, $siteBLSWithExonicFile, $outFile) = @ARGV;

#$seqRefGenomeBedFile : the bed file of the reference genome, used to extract maf alignment
#siteBLSWithExonicFile: the output of BLS pipeline, the first column is the site id, the last columns are bls scores and an indicator whether the site is overlap with an exon


print "loading regions from $seqRefGenomeBedFile ...\n" if $verbose;

my $seqs = readBedFile ($seqRefGenomeBedFile, $verbose);

my %seqHash;

foreach my $seq (@$seqs)
{
	my $seqId = $seq->{"name"};
	next if exists $seqToRemove{$seqId};

	$seqHash{$seqId} = {len=> $seq->{"chromEnd"} - $seq->{"chromStart"} + 1};
	Carp::croak "the length of sequence is too short:", Dumper ($seq), "\n" unless $seqHash{$seqId}->{"len"} > 2 * $ext;
}

my $n = keys %seqHash;

print "$n sequences loaded\n" if $verbose;



print "loading rna secondary structure data from $seqRefProbUnpairFile ...\n" if $verbose;


my $unpairWordLen = 0;

my $fin;
open ($fin, "<$seqRefProbUnpairFile") || Carp::croak "can not open file $seqRefProbUnpairFile to read\n";

while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;
	my @cols = split (/\t/, $line);
	my $seqId = shift @cols;
	
	next if exists $seqToRemove {$seqId};
	next unless exists $seqHash{$seqId};

	#if ($unpairWordLen == 0)
	#{
	$unpairWordLen = $seqHash{$seqId}->{"len"} - @cols + 1 if $unpairWordLen == 0;
	#	Carp::croak "data inconsistency in $seqRefProbUnpairFile: wordLen = $unpairWordLen, seqId=$seqId, seqLen=", $seqHash{$seqId}->{"len"}, ", pu number=", ($#cols + 1), "\n";
		
	#}
	Carp::croak "data inconsistency in $seqRefProbUnpairFile: wordLen = $unpairWordLen, seqId=$seqId, seqLen=", $seqHash{$seqId}->{"len"}, ", pu number=", ($#cols + 1), "\n" 
	if $unpairWordLen != $seqHash{$seqId}->{"len"} - @cols + 1;
	
	$seqHash{$seqId}->{"pu"} = \@cols;
	
	#Carp::croak "seqLen=", $seqHash{$seqId}->{"len"}, ", wordlen= $unpairWordLen\n";
}
close ($fin);


print "pu len = $unpairWordLen\n" if $verbose;

print "formatting site info from $siteBLSWithExonicFile ...\n" if $verbose;

my $fout;

open ($fin, "<$siteBLSWithExonicFile") || Carp::croak "cannot open file $siteBLSWithExonicFile to read\n";
open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

my $i =0;




while (my $line= <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;

	print "$i ...\n" if $i % 5000 == 0 && $verbose;
	$i++;

	my @cols = split (/\t/, $line);
	
	my $siteId = shift @cols;
	my $exonic = pop @cols;
	my $bls = pop @cols;

	my ($seqId, $mafcoord, $strand, $siteStart, $siteEnd) = split ("//", $siteId);
	$siteEnd -= 1;

	next unless exists $seqHash{$seqId};

	my $seqLen = $seqHash{$seqId}->{"len"};


	next if $siteEnd < $endToTruncate;
	next if $siteStart >= $seqLen - $endToTruncate;

	my $uifEnd = $ext - 1;
	my $exonStart = $ext;
	my $exonEnd = $seqLen - $ext - 1;
	my $difStart = $seqLen - $ext;

	my $siteLen = $siteEnd - $siteStart + 1;

	#determine how many values we need to get to average
	my $probUnpairStart = $siteStart;
	$probUnpairStart -= int (($unpairWordLen - $siteLen) / 2) if $unpairWordLen > $siteLen;
	#$probUnpairStart = 0 if $probUnpairStart < 0; 

	#Carp::croak "siteStart = $siteStart, puStart = $probUnpairStart\n";

	my $n = 1;
	$n = $siteLen - $unpairWordLen + 1 if $siteLen > $unpairWordLen;

	my $probUnpairEnd = $probUnpairStart + $n - 1;

	$probUnpairStart = 0 if $probUnpairStart < 0;
	$probUnpairEnd = 0 if $probUnpairEnd < 0;

	$probUnpairEnd = @{$seqHash{$seqId}->{"pu"}} - 1 if $probUnpairEnd >= @{$seqHash{$seqId}->{"pu"}};
	$probUnpairStart = $probUnpairEnd if $probUnpairStart > $probUnpairEnd;#this may happen when $siteLen < $unpairWordLen at the end of sequence

	#Carp::croak "i=$i, seqid = $seqId, seqlen=$seqLen, siteStart = $siteStart, siteEnd = $siteEnd, start = $probUnpairStart, end = $probUnpairEnd\n" if $probUnpairEnd < $probUnpairStart;

	my @probUnpairSite = @{$seqHash{$seqId}->{"pu"}}[$probUnpairStart..$probUnpairEnd];
	
	#Carp::croak "pu = ", join ("\t", @probUnpairSite), "\n";
	
	my $hasNa =0;
	foreach my $p (@probUnpairSite)
	{
		$hasNa = 1 if $p eq 'NA';
	}
	next if $hasNa;

	my $probUnpair = 'NA';
	$probUnpair = mean (\@probUnpairSite) unless $hasNa; 
	
	my $region = "unknown";
	if (not $keepCoord)
	{

		if ($siteStart < $exonStart)
		{
			$region = "uif";
		}
		elsif ($siteStart >= $difStart)
		{
			$region = "dif";
			$siteStart -= $difStart;
			$siteEnd -= $difStart;
		}
		else
		{
			$region = "exon";
			$siteStart -= $exonStart;
			$siteEnd -= $exonStart;
		}
	}
	$siteEnd += 1;

	print $fout join ("\t", $seqId, $siteStart, $siteEnd, $siteId, $bls, '+', $region, $bls, $probUnpair, $exonic), "\n";
}
close ($fin);
close ($fout);


