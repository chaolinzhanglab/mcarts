#!/usr/bin/perl -w
#


use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;
use Cwd;

use Bio::SeqIO;
use Data::Dumper;

use Common;
use Bed;

my $prog = basename($0);


my $verbose = 0;
my $species = "";
my $trace = 0;
my $refSpecies = 'hg18';
my $noMergeRef = 0;

my @speciesList = qw(hg18 panTro1 rheMac2 mm8 rn4 oryCun1 canFam2 bosTau2 dasNov1 loxAfr1 echTel1 monDom4 galGal2 xenTro1 tetNig1 fr1 danRer3);

GetOptions (
			'v|verbose'=>\$verbose,
			'trace'=>\$trace,
			'species:s'=>\$species,
			'ref:s'=>\$refSpecies,
			'no-merge-ref'=>\$noMergeRef,
);



#if genomeDir is specified, use that
#

if (@ARGV != 2)
{
	print "cluster motif sites detected in maf sequences\n";
	print "Usage: $prog [options] <in.bed> <out.tab>\n";
	
	print "OPTIONS:\n";
	print " -species    [string]: list of species in the maf files ($species)\n";
	print " -ref        [string]: refernece species ($refSpecies)\n";
	print " -no-merge-ref       : do not merge sites in the reference genome\n";
	print " -trace              : trace original coordinates in the reference genome if possible\n";
	print " -v                  : verbose\n";
	print "\n";
	exit (1);
}


my ($in, $out) = @ARGV;


#load species
if ($species ne '')
{
	Carp::croak "Species file $species does not exist\n" unless -f $species;
	my $fin;

	@speciesList = ();
	open ($fin, "<$species") || Carp::croak "can not open file $species to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		push @speciesList, $line;
	}
	close ($fin);
	my $ns = @speciesList;

	print "$ns species loaded\n" if $verbose;
}

my %speciesHash = map {$_=>1} @speciesList;

=obsolete
foreach my $sp (@speciesList)
{
	$speciesHash{$sp} = 1;
}
=cut

Carp::croak "The reference species $refSpecies does not exist in the list of species\n" unless exists $speciesHash{$refSpecies};

print "loading regions from $in ...\n" if $verbose;
my $regions = readBedFile($in, $verbose);

my $n = @$regions;
print "$n regions loaded\n" if $verbose;

my %siteHash;
foreach my $r (@$regions)
{
	my $chrom = $r->{"chrom"};
	my $chromStart = $r->{"chromStart"};
	my $chromEnd = $r->{"chromEnd"};
	
	my @cols = split (/\./, $chrom);
	my $sp = pop @cols;
	my $seqId = join (".", @cols);
	#my ($seqId, $sp) = split (/\./, $chrom);
	#Carp::croak "no name information in ", Dumper ($r), "\n" unless exists $r->{"name"};
	#my $name = uc ($r->{"name"}); 
	
	#Carp::croak "$chrom: $sp does not exist\n" unless exists $speciesHash{$sp};
	next unless exists $speciesHash{$sp};

	Carp::croak "no strand information in ", Dumper ($r), "\n" unless exists $r->{"strand"};
	my $strand = $r->{"strand"};

	$r->{"assembly"} = $sp;
	
	push @{$siteHash{$seqId}->{$strand}}, $r;
}


print "clustering motif sites ...\n" if $verbose;

my $fout;
open ($fout, ">$out") || Carp::croak "can not open file $out to write\n";

print $fout join ("\t", "SEQ", @speciesList), "\n";

my $total = 0;
foreach my $seqId (sort keys %siteHash)
{
	#print "$seqId\n";
	my $siteInSeq = $siteHash{$seqId};
	foreach my $strand (qw(+ -))
	{
		next unless exists $siteInSeq->{$strand};
		
		my $siteInSeqStrand = $siteInSeq->{$strand};

		#print Dumper ($siteInSeqStrand), "\n";
		
		#my $nsite = @$siteInSeqStrand;
	
		#building connect matrix directly require too much memory
		#so we cluster overlapping sites first
		#the we refine and identify orthologous sites

        #definition: clusterRegions ($regionsOnChrom, $strand, $maxGap, $overlapFraction, $collapse); 

		my $overlapSiteClusters = clusterRegions ($siteInSeqStrand, $strand, 0, 0, 0); 
		#changed on 12/18/2010 due to changes in the subroutine
		
		foreach my $clust (@$overlapSiteClusters)
		{
			#define orthologous sites (more strict definition)in a cluster
			my @siteInClust = @$siteInSeqStrand[@$clust];
			
			my $connect = [];
			my $nsite = @siteInClust;
		
			#$total += $nsite;
			for (my $i = 0; $i < $nsite; $i++)
			{
				my $r1 = $siteInClust[$i];
				$connect->[$i][$i] = 0;
			
				for (my $j = $i+1; $j < $nsite; $j++)
				{
					my $r2 = $siteInClust[$j];
					#print "r1=", Dumper($r1), "\nr2=", Dumper ($r2), "\n";
					if (($r1->{"chromStart"} <= $r2->{"chromStart"} && $r1->{"chromEnd"} >= $r2->{"chromEnd"})
					||  ($r1->{"chromStart"} >= $r2->{"chromStart"} && $r1->{"chromEnd"} <= $r2->{"chromEnd"}))
					{
						$connect->[$i][$j] = 1;
					}
					else
					{
						$connect->[$i][$j] = 0;
					}
					$connect->[$j][$i] = $connect->[$i][$j];
				
				}
			}

			#for (my $i = 0; $i < $nsite; $i++)
			#{
			#	print join (" ", @{$connect->[$i]}), "\n";
			#}
			my $orthologousSites = Common::matrix2clusters ($connect);

			foreach my $s (@$orthologousSites)
			{
				my $n = @$s; #how many species have the site

				#print join ("\t", @$c), "\n";
				my @regions = @siteInClust[@$s]; #the cordinates and other information
				$total += @regions;	
				#Carp::croak Dumper (\@regions), "\n";
				printOrthologousSites ($fout, \@regions, $seqId, $strand);
			} #end of a site
		} #end of a cluster
	}#end of  a strand
}#end of a sequence


print "total site = $total\n";
close ($fout);

#used global variables, should never be moved to packages without modification
sub printOrthologousSites
{
	my ($fout, $regions, $seqId, $strand) = @_;

	my @regions = @$regions;

	my $chromStart = $regions[0]->{"chromStart"};
	my $chromEnd = $regions[0]->{"chromEnd"};
	my %s;	#number of sites present a particular genome

	my @siteInRef;
			
	my $traceInfo = "";
			
	foreach my $r (@regions)
	{
		$chromStart = $r->{"chromStart"} if $chromStart > $r->{"chromStart"};
		$chromEnd = $r->{"chromEnd"} if $chromEnd < $r->{"chromEnd"};
		if (exists $s{$r->{"assembly"}})
	   	{
			$s{$r->{"assembly"}}++;
		}
		else
		{
			$s{$r->{"assembly"}} = 1;
		}

		if ($r->{"assembly"} eq $refSpecies)
		{
			$traceInfo .= "//" . $r->{"name"} if $trace;
			push @siteInRef, $r;
		}
		#$traceInfo .= "//" . $r->{"name"} if ($trace && $r->{"assembly"} eq $refSpecies);
	}
	
	my @counts;

	$s{$refSpecies} = 1 if $noMergeRef && @siteInRef > 0;
			
	foreach my $assembly (@speciesList)
	{
		my $ct = 0;
		$ct = $s{$assembly} 
		if (exists $s{$assembly});
		push @counts, $ct;
	}
			
	if ($noMergeRef)
	{
		if (@siteInRef > 0)
		{
			foreach my $r (@siteInRef)
			{
				my $chromStart = $r->{"chromStart"};
				my $chromEnd = $r->{"chromEnd"} + 1;
				my $name = "$seqId//$chromStart-$chromEnd//$strand";
				$name .= "//" . $r->{"name"} if $trace;
				print $fout join ("\t", $name, @counts), "\n";
			}
		}
		else
		{
			$chromEnd += 1;
			my $name = "$seqId//$chromStart-$chromEnd//$strand";
			print $fout join ("\t", $name, @counts), "\n";
		}
	}
	else
	{
		$chromEnd += 1;
		my $name = "$seqId//$chromStart-$chromEnd//$strand";
		$name .= $traceInfo if $trace;
		print $fout join ("\t", $name, @counts), "\n";
	}
}
