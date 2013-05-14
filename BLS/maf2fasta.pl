#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;

use Sequence;


my $prog = basename($0);

my $noGap = 0;
my $noMap = 0;
my $maskRepeats = 0;
my $verbose = 0;
my $segment = 0;

GetOptions ('m|mask_repeat'=>\$maskRepeats,
			'G|no-gap'=>\$noGap,
			'nomap'=>\$noMap,
			's|segment-sequence:i'=>\$segment,
			'v|verbose'=>\$verbose
);

if (@ARGV != 2)
{
	print "Convert MAF files to a fasta file\n";
	print "Usage: $prog [options] <in.maf> <out.fa>\n";
	print "          use \"-\" for STDIN or STDOUT\n";
	print "          gzip compressed input is allowed\n";
	print "OPTIONS:\n";
	print " -m      : mask repeats <on|[off]>\n";
	print " -G      : remove gap <on|[off]>\n";
	print " -nomap  : do not write gap map\n";
	print " -s [int]: segment sequences ($segment) (0=no segment)\n";
	print " -v      : verbose <on|[off]>\n";
	print "\n";
	exit (1);
}

my ($inMafFile, $outFastaFile) = @ARGV;


my $fin;

if ( $inMafFile eq "-")
{
    $fin = *STDIN;
}
else
{
	if ($inMafFile =~/\.gz$/)
	{
		open ($fin, "gunzip -c $inMafFile | ") || Carp::croak "cannot open file $inMafFile to read\n";
	}
	else
	{
		open ($fin, "<$inMafFile") || Carp::croak "can not open file $inMafFile to read\n";
	}
}
my $fout;

if ($outFastaFile eq '-')
{
	$fout = *STDOUT;
}
else
{
	open ($fout, ">$outFastaFile") || Carp::croak "can not open file $outFastaFile to write\n";
}

my $blockName = "";
my $seqPrefix = "";
my $blockIter = -1;
while (my $line = <$fin>)
{
	chomp $line;
	if ($line =~/^\#BLOCK_NAME=(.*?)$/)
	{
		$blockName = $1;
	}
	elsif ($line =~/^a/)
	{
		$blockIter++;
		if ($blockName ne "")
		{
			$seqPrefix = $blockName;
		}
		else
		{
			$seqPrefix = "b$blockIter";
		}
	}
	elsif ($line =~/^\s*$/)
	{
		$blockName = "";
		print $fout "\n";
	}
	elsif ($line =~/^s/)
	{
		my @cols = split (/\s+/, $line);
		my $db = $cols[1];
		if ($db=~/^(\w+)\./)
		{
			$db = $1;
		}
		my $seq = $cols[6];

		my $header = "$seqPrefix.$db";
		if ($noGap)
		{
			my $gapLen = 0;
			my @gapMap;
			if ($seq!~/^\-/)
			{
				push @gapMap, "0->0";
			}
			while ($seq=~/(\-+)/g)
			{
				$gapLen += length($1);

				#print "$seqPrefix.$db=$gapLen\n";
				my $gapEnd = pos($seq);
				my $gaplessEnd = $gapEnd - $gapLen;
				
				push @gapMap, "$gaplessEnd->$gapEnd";
			}
			$header .= "\tMAP=" . join (",", @gapMap) unless $noMap;
			
			$seq =~s/-//g;
		}

		if ($maskRepeats)
		{
			$seq =~tr/acgt/nnnn/;
		}
		print $fout ">$header\n";
		if ($segment > 0)
		{
			print $fout segmentStr ($seq, int($segment)), "\n";
		}
		else
		{
			print $fout $seq, "\n";
		}
	}
}

close ($fin) if $inMafFile ne '-';
close ($fout) if $outFastaFile ne '-';


