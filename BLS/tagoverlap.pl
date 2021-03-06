#!/usr/bin/perl -w


use strict;
use warnings;

use Getopt::Long;
use Carp;
use File::Basename;
use Data::Dumper;

use MyConfig;
use Bed;
use Common;

my $prog = basename ($0);

my $regionFile = "";
my $separateStrand = 0;
my $completeOverlap = 0;
my $verbose = 0;
my $ext5 = 0;
my $ext3 = 0;
my $cache = getDefaultCache ($prog);
my $delimitor = "//";
my $keepScore = 0;
my $keepCache = 0;

my $reverse = 0; #find nonoverlapping tags
my $big = 0;
my $minBlockSize = 2000000;

GetOptions (
		'ext5:i'=>\$ext5,
		'ext3:i'=>\$ext3,
		'region:s'=>\$regionFile,
		'ss|sep-strand'=>\$separateStrand,
		'complete-overlap'=>\$completeOverlap,
		'r'=>\$reverse,
		'd:s'=>\$delimitor,
		'big'=>\$big,
		'minBlockSize:i'=>\$minBlockSize,
		'keep-score'=>\$keepScore,
		'c|cache:s'=>\$cache,
		'keep-cache'=>\$keepCache,
		'v|verbose'=>\$verbose);

if (@ARGV != 2)
{
	print "Find tags overlap with particular regions\n";
	print "Usage1: $prog [options] <tags.bed> <overlap.out>\n";
	print "Usage2: $prog [options] <tags_dir_by_chrom> <overlap.out>\n";
	print "OPTIONS:\n";
	print " -big              : either region or tag file is big\n";
	print " -ext5        [int]: extension of tags at the 5' end\n";
	print " -ext3        [int]: extension of tags at the 3' end\n";
	print " -region     [file]: a bed file with regions to count tag numbers.\n";
	print " -ss               : separate strand (off)\n";
	print " --complete-overlap: requires complete overlap of the tag with the region\n";
	print " -r                : reverse mode to print tags without ovrlap with the region (off)\n";
	print " -d        [string]: delimitor to separate ids\n";
	print " --keep-score      : keep tag score\n";
	print " -c        [string]: cache dir ($cache)\n";
	print " --keep-cache      : keep cache when the job is done\n";
	print " -v                : verbose\n";
	exit (1);
}

my ($tagBedFile, $outFile) = @ARGV;

print "keep-score = $keepScore\n" if $verbose;
print "loading tags from $tagBedFile ...\n" if $verbose;
my %tagCount;

if (-f $tagBedFile)
{
	if ($big)
	{
		die "$cache already exist\n" if (-d $cache);
		die "$cache already exist\n" if (-f $cache);
		system ("mkdir $cache");
		my $tagCache = "$cache/tags";
		system ("mkdir $tagCache");

		my $ret = splitBedFileByChrom ($tagBedFile, $tagCache, v=>$verbose, "sort"=>1);
		%tagCount = %$ret;
	}
	else
	{
		my $tags = readBedFile ($tagBedFile, $verbose);

		foreach my $t (@$tags)
		{
			my $chrom = $t->{'chrom'};
			push @{$tagCount{$chrom}}, $t;
		}
	}
}
elsif (-d $tagBedFile)
{
	print "check tags in dir $tagBedFile ...\n" if $verbose;

	#we assume these files are already sorted
	my @fnames = `ls $tagBedFile`;
	chomp @fnames;

	foreach my $f (@fnames)
	{
		if ($f =~/^(.*?)\.bed$/i)
		{
			my $chrom = $1;
			#print "$tagBedFile/$f", "\n";
			my $n = `wc -l $tagBedFile/$f`;
			chomp $n;
			$n=~/^(\d+)\s/;
			$n = $1;
			$tagCount{$chrom} = {f=>"$tagBedFile/$f", n=>$n};
		}
	}
}
else
{
	Carp::croak "invalid tag file or dir: $?\n";
}

print "get tag count broken down into chromosomes ...\n" if $verbose;

foreach my $chrom (sort keys %tagCount)
{
	if (ref($tagCount{$chrom}) eq 'HASH')
	{
		my $n = $tagCount{$chrom}->{"n"};
		print "$chrom : $n tags\n" if $verbose;
	}
	else
	{
		my $n = @{$tagCount{$chrom}};
		print "$chrom : $n tags\n" if $verbose;
	}
}

my %regionHash= ();

Carp::croak "$regionFile does not exists\n" unless -f $regionFile;
print "reading regions from $regionFile ...\n" if $verbose;

my $n = 0;
if ($big)
{
	my $regionCache = "$cache/regions";
	system ("mkdir $regionCache");

	my $ret = splitBedFileByChrom ($regionFile, $regionCache, $verbose);
	%regionHash = %$ret;

	foreach my $chrom (keys %regionHash)
	{
		$n += $regionHash{$chrom}->{"n"};
	}
}
else
{
	my $regions = readBedFile ($regionFile, $verbose);
	$n = @$regions;	
	foreach my $r (@$regions)
	{
		my $chrom = $r->{'chrom'};
		push @{$regionHash{$chrom}}, $r;
	}
}

print "$n regions loaded\n" if $verbose;



print "searching for overlapping tags ...\n" if $verbose;

my @strand = ('b');
@strand = qw(+ -) if $separateStrand;

my $fout;
open ($fout, ">$outFile") || Carp::croak "can not open file $outFile to write\n";


foreach my $chrom (sort keys %regionHash)
{
	print "processing chrom $chrom ...\n" if $verbose;
	next unless exists $tagCount {$chrom};

	#prepare regions

	print "preparing regions ...\n" if $verbose;
	my $regionsOnChrom;

	#Carp::croak ref($regionHash{$chrom}), "\n";
	if (ref($regionHash{$chrom}) eq 'ARRAY')
	{
		$regionsOnChrom = $regionHash{$chrom};
	}
	else
	{
		my $tmpFile = $regionHash{$chrom}->{'f'};
		print "loading regions on chromsome $chrom from $tmpFile...\n" if $verbose;
		$regionsOnChrom = readBedFile ($tmpFile, $verbose);
		my $n = @$regionsOnChrom;
		print "$n regions loaded on chromosome $chrom\n" if $verbose;
	}

	my %regionsOnChromHash;

	foreach my $r (@$regionsOnChrom)
	{
		$r->{"score"} = 0;
		$r->{"name"} = $r->{"chrom"} . ":" . ($r->{"chromStart"} + 1) . "-" . ($r->{"chromEnd"} + 1) unless exists $r->{"name"};

		if ($separateStrand)
		{
			Carp::croak "no strand information in ", Dumper ($r), "\n" unless exists $r->{"strand"};
			my $s = $r->{"strand"};
			
			push @{$regionsOnChromHash{$s}}, $r;
		}
		else
		{
			push @{$regionsOnChromHash{'b'}}, $r;
		}
	}

	#prepare tags
	
	print "preparing tags ...\n" if $verbose;

	if (ref ($tagCount{$chrom}) eq 'ARRAY')
	{
		my $tagsOnChrom = $tagCount{$chrom};
		processBlock ($tagsOnChrom, \%regionsOnChromHash, $chrom, $separateStrand, $ext5, $ext3, $reverse, $fout);
	}
	else
	{
		my $tmpFile = $tagCount{$chrom}->{"f"};
		my $fin;

		open ($fin, "<$tmpFile") || Carp::croak "cannot open file $tmpFile to read\n";
		my $iter = 0;
		while (my $block = readNextBedBlock ($fin, max (1, $ext5 + $ext3 + 1), 0, minBlockSize=> $minBlockSize))
		{
			my $n = @$block;
			print "block $iter : $n tags\n" if $verbose;
			$iter++;
			processBlock ($block, \%regionsOnChromHash, $chrom, $separateStrand, $ext5, $ext3, $reverse, $fout);
		}
		close ($fin);
	}
}
close ($fout);

system ("rm -rf $cache") unless $keepCache;


sub processBlock
{
	my ($tagsOnChrom, $regionsOnChromHash, $chrom, $separateStrand, $ext5, $ext3, $reverse, $fout) = @_;
	
	my @strand = ('b');
	@strand = qw(+ -) if $separateStrand;

	my %tagsOnChromHash;

	foreach my $t (@$tagsOnChrom)
	{

		if ($separateStrand)
		{
			Carp::croak "no strand information in ", Dumper ($t), "\n" unless exists $t->{"strand"};
			my $s = $t->{"strand"};
			
			push @{$tagsOnChromHash{$s}}, $t;
		}
		else
		{
			push @{$tagsOnChromHash{'b'}}, $t;
		}
		
		$t->{"strand"} = '+' unless exists $t->{"strand"};
	
		#extend tags if necessary
		if ($t->{"strand"} eq '+')
		{
			$t->{"chromStart"} -= $ext5;
			$t->{"chromEnd"} += $ext3;
		}
		else
		{
			$t->{"chromStart"} -= $ext3;
			$t->{"chromEnd"} += $ext5;
		}
	}


	foreach my $s (@strand)
	{
		print "finding overlapping tags on strand $s of $chrom ...\n" if $verbose;
		next unless exists $regionsOnChromHash->{$s} && exists $tagsOnChromHash{$s};

		my $nr = @{$regionsOnChromHash->{$s}};
		my $nt = @{$tagsOnChromHash{$s}};

		print "region number = $nr, tag number = $nt\n" if $verbose;

		findOverlapTags ($tagsOnChromHash{$s}, $regionsOnChromHash->{$s}, $reverse, $fout);
		
		if ($reverse)
		{
			foreach my $t (@{$tagsOnChromHash{$s}})
			{
				print $fout bedToLine ($t), "\n" unless exists $t->{"overlap"} && $t->{"overlap"} == 1;
			}
		}
	}
}

sub findOverlapTags
{
	my ($tagsUnSorted, $regionsUnSorted, $reverse, $fout) = @_;

	my @tags = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$tagsUnSorted;
	my $tags = \@tags;

	my $chromEndExt = 0;
	foreach my $t (@$tags)
	{
		$chromEndExt = $t->{"chromEnd"} if $t->{"chromEnd"} > $chromEndExt;
		$t->{"chromEndExt"} = $chromEndExt;
	}
	
	my @regions = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$regionsUnSorted;
	my $regions = \@regions;
	
	my $firstTagIdx = 0; #the first tag that overlap with or on the right of the current window

	foreach my $r (@$regions)
	{
		
		my $chromStart = $r->{"chromStart"};
		my $chromEnd = $r->{"chromEnd"};
		
		while ($firstTagIdx < @$tags && $tags->[$firstTagIdx]->{"chromEndExt"} < $chromStart)
		{
			$firstTagIdx++;
		}
		
		my $i = $firstTagIdx;
		
		while ($i < @$tags && $tags->[$i]->{"chromStart"} <= $chromEnd)
		{
			if ($tags->[$i]->{"chromEnd"} >= $chromStart)
			{
				my $overlap = 1;
				if (exists $tags->[$i]->{"blockCount"})
				{
					#check if the region is in the intron of the tags
					for (my $k = 1; $k < $tags->[$i]->{"blockCount"}; $k++)
					{
						my $intronStart = $tags->[$i]->{"chromStart"} + $tags->[$i]->{"blockStarts"}->[$k-1] + $tags->[$i]->{"blockSizes"}->[$k-1];
						my $intronEnd = $tags->[$i]->{"chromStart"} + $tags->[$i]->{"blockStarts"}->[$k] - 1;

						if ($intronStart <= $chromStart && $intronEnd >= $chromEnd)
						{
							$overlap = 0;
							last;
						}
					}
				}
				if ($overlap)
				{
					$tags->[$i]->{"overlap"} = 1;
					if ($reverse == 0)
					{
						my %r = %{$tags->[$i]};
						$r{"name"} .= $delimitor . $r->{"name"};
						my $start = max($chromStart, $tags->[$i]->{"chromStart"});
						my $end = min ($chromEnd, $tags->[$i]->{"chromEnd"});
	
						my $score = ($end - $start + 1) / ($tags->[$i]->{"chromEnd"} - $tags->[$i]->{"chromStart"} + 1);

						#$tags->[$i]->{"name"} .= "//" . $r->{"name"};
						$r{"score"} = $score unless $keepScore;
						
						if ($completeOverlap)
						{
							print $fout bedToLine (\%r), "\n" if $score == 1;
						}
						else
						{
							print $fout bedToLine (\%r), "\n";
						}
					}
				}
			}
			$i++;
		}
	}
}


