#!/usr/bin/perl -w
#


use lib "/mnt/raygun/zhangc/perl_lib2";	#own library
use lib "/mnt/raygun/zhangc/perl_lib"; #for Bioperl

use strict;
use File::Basename;
use Getopt::Long;
use Carp;
use Cwd;
use Bio::SeqIO;
use Data::Dumper;

use Bed;


my $prog = basename($0);

my $trace = 0;
my $debug = 0;
my $verbose = 0;

GetOptions (
			'd|debug'=>\$debug,
			'trace'=>\$trace,
			'v|verbose'=>\$verbose
);



#if genomeDir is specified, use that
#

if (@ARGV != 3)
{
	print "Map gapless coordinates back to coordinates with gaps in the MAF file\n";
	print "Usage: $prog [options] <in.bed> <in.maf2fa> <out.bed>\n";
	print " <in.bed>    : the coordinates according to gapless sequences from the MAF file\n";
	print " <in.maf2fa> : the gapless sequences extracted from the MAF file (gzip file is allowed)\n";
	print " <out.bed>   : the coordinates counting gaps in the map file\n";
	print "OPTIONS:\n";
	print " -debug      : print debug information\n";
	print " -trace      : trace original coordinates (written in the name column)\n";
	print " -v          : verbose\n";
	print "\n";
	exit (1);
}


my ($inBedFile, $inFastaFile, $outBedFile) = @ARGV;


my %gapMap;

print "loading gap map information from $inFastaFile ...\n" if $verbose;
my $fin;

if ($inFastaFile =~/\.gz$/)
{
	open ($fin, "gunzip -c $inFastaFile | ") || Carp::croak "cannot open file $inFastaFile to read\n";
}
else
{
	open ($fin, "<$inFastaFile") || Carp::croak "can not open file $inFastaFile to read\n";
}
while (my $line = <$fin>)
{
	chomp $line;
	next unless $line=~/\>/;
	my ($seqName, $map) = split (/\s+/, $line);
	$seqName =~/^\>(.*?)$/;
	$seqName = $1;
	if ($map=~/^MAP=(.*?)$/)
	{
		$map = $1;
	}
	else
	{
		$map = "";
	}
	my @entries = split (/\,/, $map);
	foreach my $e (@entries)
	{
		my ($from, $to) = split (/\-\>/, $e);
		push @{$gapMap{$seqName}}, {from=>$from, to=>$to};
	}
}
close ($fin);

my $nmap = keys %gapMap;
print "$nmap sequences have gap maps loaded\n" if $verbose;



print "loading regions from $inBedFile ...\n";
my $regions = readBedFile($inBedFile, $verbose);


my $n = @$regions;

my %regionHash;

#sort out according to gene contigs
foreach my $r (@$regions)
{
	my $chrom = $r->{"chrom"};
	push @{$regionHash{$chrom}}, $r;	
}


print "$n regions loaded\n" if $verbose;


print "mapping coordinates ...\n" if $verbose;

#my @newRegions;

my $fout;

open ($fout, ">$outBedFile") || Carp::croak "can not open file $outBedFile to write\n";

my $i = 0;

foreach my $chrom (sort keys %regionHash)
{
	if (not exists $gapMap {$chrom})
	{
		warn "No gap map for seq $chrom\n";
		next;
	}
	my $regionsOnChrom = $regionHash{$chrom};
	my @regionsOnChromSort = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$regionsOnChrom;

	my $map = $gapMap{$chrom};
	print Dumper ($map), "\n" if $debug;

	print "processing $chrom ...\n" if $verbose || $debug;
	my $mapStartIdx = 0;

	foreach my $r (@regionsOnChromSort)
	{
		print "$i of $n ...\n" if $i % 1000 == 0 && ($verbose || $debug);
		$i++;

		#my $chrom = $r->{"chrom"};
		my $chromStart = $r->{"chromStart"};
		my $chromEnd = $r->{"chromEnd"};
	
		if ($trace)
		{
			$r->{"name"} = "$chromStart//" . ($chromEnd + 1);
		}

		my $map = $gapMap{$chrom};

		for (; $mapStartIdx < @$map; $mapStartIdx++)
		{
			my $m = $map->[$mapStartIdx];
			last if $m->{'from'} > $chromStart;
		}

		my $j = $mapStartIdx;
		if ($mapStartIdx < @$map)
		{
			print "j =$j, chromStart = $chromStart, last map entry=", $map->[$j-1]->{"from"}, ", current map entry=", $map->[$j]->{"from"}, "\n" if $debug;
		}
		else
		{
			print "j = $j, chromStart = $chromStart, last map entry=", $map->[$j-1]->{"from"}, ", last block\n" if $debug;
		}

		my $mStart = $map->[$j-1];
		Carp::croak "chromStart = $chromStart, blockStart = ", $mStart->{'from'}, ", inconsistent\n" unless $chromStart >= $mStart->{'from'};
   	
		my $chromStartNew = $mStart->{"to"} + $chromStart - $mStart->{'from'};
	
		for (; $j < @$map; $j++)
		{
			my $m = $map->[$j];
			last if $m->{'from'} > $chromEnd;
		}
	
		if ($j < @$map)
		{
			print "j = $j, chromEnd = $chromEnd, last map entry=", $map->[$j-1]->{"from"}, ", current map entry=", $map->[$j]->{"from"}, "\n" if $debug;
		}
		else
		{
			print "j = $j, chromEnd = $chromEnd, last map entry=", $map->[$j-1]->{"from"}, ", last block\n" if $debug;
		}
		my $mEnd = $map->[$j-1];
		Carp::croak "chromStart = $chromEnd, blockStart = ", $mEnd->{'from'}, ", inconsistent\n" unless $chromEnd >= $mEnd->{'from'};
	
		my $chromEndNew = $mEnd->{"to"} + $chromEnd - $mEnd->{'from'};

		$mapStartIdx--;

		$r->{"chromStart"} = $chromStartNew;
		$r->{"chromEnd"} = $chromEndNew;
	
		print $fout bedToLine ($r), "\n";
		#push @newRegions, $r;
	}
	@regionsOnChromSort = ();
}
close ($fout);

