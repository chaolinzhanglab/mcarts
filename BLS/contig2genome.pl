#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Data::Dumper;

use Bed;
use Wiggle;

my $prog = basename ($0);

my $verbose = 0;
my $format = 'bed';
my $trackName = '';
my $uniqId = 0;
my @ARGV0 = @ARGV;

GetOptions ('v|verbose'=>\$verbose,
		'uniq-id'=>\$uniqId,
		'f|format:s'=>\$format,
		'n|name:s'=>\$trackName);

if ($format eq 'wig')
{
	Carp::croak "Track name must be specified for wig format\n" if $trackName eq '';
}

if (@ARGV != 3)
{
	print "convert contig coordiantes to genome coordinates\n";
	print "Usage: $prog [options] <contig.bed> <ts_on_contig.bed> <out.bed>\n";
	print "OPTIONS\n";
	print " -v  : verbose (off)\n";
	print " --uniq-id  : make uniq transcript id\n";
	print " -f [string]: input format ([bed]|wig)\n";
	print " -n [string]: track name\n";
	exit (1);
}

my ($contigBedFile, $tsOnContigFile, $tsOnGenomeFile) = @ARGV;


print "CMD=$prog ", join(" ", @ARGV0), "\n" if $verbose;


print "reading contig bed file from $contigBedFile ...\n" if $verbose;
my $contigs = readBedFile ($contigBedFile);
my $ncontig = @$contigs;

print "$ncontig loaded\n" if $verbose;

my %contigHash;
foreach my $c (@$contigs)
{
	$contigHash{$c->{"name"}} = $c;
	$c->{"strand"} = '+' unless exists $c->{"strand"};
}

if ($format eq 'bed')
{
	print "reading bed file from $tsOnContigFile ...\n" if $verbose;
	my $tsOnContig = readBedFile ($tsOnContigFile, $verbose);
	my $nts = @$tsOnContig;

	print "$nts rows loaded\n" if $verbose;

	#print Dumper ($tsOnContig), "\n";

	my @tsOnGenome;

	print "coordinate conversion ...\n" if $verbose;

	my $i = 0;
	foreach my $ts (@$tsOnContig)
	{
		if ($verbose)
		{
			print "processing rows $i...\n" if $i % 5000 == 0;
		}
		my $tsId = $ts->{"name"};
		$tsId = $i . "#" . $tsId if $uniqId;
		
		$i++;
		my $contigId = $ts->{"chrom"};

		if (not exists $contigHash{$contigId})
		{
			warn "$tsId: contig $contigId does not exists\n";
			next;
		}
		my $contig = $contigHash{$contigId};

		#print Dumper ($contig), "\n";
		my %tsCpy = ();
		%tsCpy = %$ts;

		#make it a complete copy, important! took me two hours to fix this
		if (exists $ts->{"blockCount"})
		{
			my @blockSizes = @{$ts->{"blockSizes"}};
			my @blockStarts = @{$ts->{"blockStarts"}};
		
			$tsCpy{"blockSizes"} = \@blockSizes;
			$tsCpy{"blockStarts"} = \@blockStarts;
		}
	
	
	
		my $tsStartOnContig = $ts->{"chromStart"};
		my $tsEndOnContig = $ts->{"chromEnd"};

	
		my $tmp = contigToGenome ($contig, $tsStartOnContig, $tsEndOnContig);
	
		$tsCpy{'chrom'} = $contig->{"chrom"};
		$tsCpy{'chromStart'} = $tmp->{"chromStart"};
		$tsCpy{'chromEnd'} = $tmp->{"chromEnd"};
		$tsCpy{'strand'} = $contig->{"strand"};
		if (exists $ts->{"strand"})
		{
			$tsCpy{'strand'} = ($contig->{"strand"} eq $ts->{"strand"})? '+' : '-';
		}

		if (exists $ts->{"thickStart"})
		{
			$tsCpy{'thickStart'} = $tsCpy{'chromStart'};
			$tsCpy{'thickEnd'} = $tsCpy{'chromEnd'};
		}
	
	
		if (exists $ts->{"blockCount"})
		{
			my $blockCount = $ts->{"blockCount"};
			for (my $i = 0; $i < $blockCount; $i++)
			{
				my $blockStartOnContig = $ts->{"blockStarts"}->[$i];
				my $blockEndOnContig = $ts->{"blockStarts"}->[$i] + $ts->{"blockSizes"}->[$i] - 1;
			
				$blockStartOnContig += $tsStartOnContig;
				$blockEndOnContig += $tsStartOnContig;

			
				my $tmp = contigToGenome ($contig, $blockStartOnContig, $blockEndOnContig);
			

				if ($contig->{"strand"} eq '+')
				{
					$tsCpy{'blockStarts'}->[$i] = $tmp->{"chromStart"} - $tsCpy{'chromStart'};
				}
				else
				{
					$tsCpy{'blockStarts'}->[$blockCount - 1 - $i] = $tmp->{"chromStart"} - $tsCpy{'chromStart'};
					$tsCpy{'blockSizes'}->[$blockCount - 1 - $i] = $ts->{"blockSizes"}->[$i];
				}
			}
		}
		#print Dumper (\%tsCpy), "\n";
		push @tsOnGenome, \%tsCpy;
	}
	print "dumping output...\n" if $verbose;
	my $header = "";
	if ($trackName ne '')
	{
		$header = "track name=\"$trackName\"";
	}
	writeBedFile (\@tsOnGenome, $tsOnGenomeFile, $header);
}
elsif ($format eq 'wig')
{
	print "reading wig file from $tsOnContigFile ...\n" if $verbose;
	my $tsOnContig = readBedGraphFile ($tsOnContigFile);
	my $nts = @$tsOnContig;

	print "$nts rows loaded\n" if $verbose;

	#print Dumper ($tsOnContig), "\n";

	my %tsOnGenome;
	#my @tsOnGenome;
	
	print "coordinate conversion ...\n" if $verbose;

	my $i = 0;
	foreach my $ts (@$tsOnContig)
	{
		if ($verbose)
		{
			print "processing row $i...\n" if ($i - int($i/500)*500 == 0);
		}
		$i++;
		my $contigId = $ts->{"chrom"};
		if (not exists $contigHash{$contigId})
		{
			warn "contig $contigId does not exists\n";
			next;
		}
	
		my $contig = $contigHash{$contigId};

		#print Dumper ($contig), "\n";
		my %tsCpy = ();
		%tsCpy = %$ts;

		my $tsStartOnContig = $ts->{"chromStart"};
		my $tsEndOnContig = $ts->{"chromEnd"};

		my $tmp = contigToGenome ($contig, $tsStartOnContig, $tsEndOnContig);
		
		$tsCpy{'chrom'} = $contig->{"chrom"};
		$tsCpy{'chromStart'} = $tmp->{"chromStart"};
		$tsCpy{'chromEnd'} = $tmp->{"chromEnd"};
		my $chrom = $tsCpy{'chrom'};
		push @{$tsOnGenome{$chrom}->{$contigId}->{"blocks"}}, \%tsCpy;
		
		#push @{$tsOnGenome{$chrom}}, \%tsCpy;
		
		#push @tsOnGenome, \%tsCpy;
	}
	
	print "split data into compatible tracks ...\n" if $verbose;
	my @contigIdsAllRounds;
	
	foreach my $chrom (sort keys %tsOnGenome)
	{
		print "processing chromosome $chrom ...";
		my $tsOnChrom = $tsOnGenome{$chrom};
		#record the start and end of each gene
		foreach my $contigId (keys %$tsOnChrom)
		{
			my $blocks = $tsOnChrom->{$contigId}->{"blocks"};
			#sort blocks
			my @blocks = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$blocks;
			$tsOnChrom->{$contigId}->{"blocks"} = \@blocks;

			#record start and end point of each contig
			my $blockCount = @blocks;
			$tsOnChrom->{$contigId}->{"chromStart"} = $blocks[0]->{"chromStart"};
			$tsOnChrom->{$contigId}->{"chromEnd"} = $blocks[$blockCount - 1]->{"chromEnd"};
		}

		#sort all contigs on the chromosome
		my @contigIds = sort {$tsOnChrom->{$a}->{"chromStart"} <=> $tsOnChrom->{$b}->{"chromStart"}} keys %$tsOnChrom;

		my $roundId = 0;
		while (@contigIds > 0)
		{
			my @contigIdsThisRound;
			my @contigIdsNextRound;

			my $ncontigs = @contigIds;

			my $prevContigId = $contigIds[0];
			push @contigIdsThisRound, $contigIds[0];
			
			for (my $i = 1; $i < $ncontigs; $i++)
			{
				my $prevContig = $tsOnChrom->{$prevContigId};
				my $currContig = $tsOnChrom->{$contigIds[$i]};
				
				if ($currContig->{"chromStart"} > $prevContig->{"chromEnd"})#compatible
				{
					push @contigIdsThisRound, $contigIds[$i];
					$prevContigId = $contigIds[$i];
				}
				else#not compatible
				{
					push @contigIdsNextRound, $contigIds[$i];
				}
			}
			
			push @{$contigIdsAllRounds[$roundId]->{$chrom}}, @contigIdsThisRound;
			$roundId++;
			@contigIds = @contigIdsNextRound;
		}
		print "split into $roundId tracks\n" if $verbose;
	}
	
	my $nRound = @contigIdsAllRounds;
	print "$nRound tracks are required to make the data compatible with wig format\n" if $verbose;
	
	print "dumping output ...\n" if $verbose;
	#my $fout = new FileHandle;
	#open ($fout, ">$tsOnGenomeFile") || Carp::croak "can not open file $tsOnGenomeFile to write\n";

	my $foutId;
	open ($foutId, ">$tsOnGenomeFile.bed") || Carp::croak "can not open file $tsOnGenomeFile.bed to write\n";
	for (my $i = 0; $i < $nRound; $i++)
	{
		my $fout;
		open ($fout, ">$tsOnGenomeFile.$i") || Carp::croak "can not open file $tsOnGenomeFile.$i to write\n";
		my $contigIdsThisRound = $contigIdsAllRounds[$i];
		my $trackNameThisRound = "$trackName" . "_$i";
		my $header = "track type=wiggle_0 name=\"$trackNameThisRound\" maxHeightPixels=50:50:10";
		printContigsInWig (\%tsOnGenome, $contigIdsThisRound, $header, $fout, $foutId);

		close ($fout);
	}
	#close ($fout);
	close ($foutId);
}
else
{
	Carp::croak "incorrect format $format\n";
}



sub printContigsInWig
{
	my ($tsOnGenome, $contigIds, $header, $fout, $foutId) = @_;
	#print $fout $header, "\n";

	foreach my $chrom (sort keys %$contigIds)
	{
		my $contigIdsOnChrom = $contigIds->{$chrom};
		foreach my $contigId (@$contigIdsOnChrom)
		{
			my $contig = $tsOnGenome->{$chrom}->{$contigId};
			my $blocks = $contig->{"blocks"};
			foreach my $b (@$blocks)
			{
				print $fout join ("\t", $b->{"chrom"},
						$b->{"chromStart"},
						$b->{"chromEnd"} + 1,
						$b->{"score"}), "\n";
				print $foutId join ("\t", $b->{"chrom"}, 
						$b->{"chromStart"}, 
						$b->{"chromEnd"} + 1,
						$contigId,
						$b->{"score"}), "\n";
						
			}
		}
	}
}




