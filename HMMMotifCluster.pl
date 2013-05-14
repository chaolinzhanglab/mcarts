#!/usr/bin/perl -w

#------------------------------------------------------------------------------
# Clustered accessible RBP target sites (carts)
#
# predict RBP motif clusters by integrating site spacing, sequence conservation and accessibility
# Chaolin Zhang (czhang@rockefeller.edu), All right reserved
# created: 10/3/2008
# last update: 05/24/2011
#
#------------------------------------------------------------------------------
#


use strict;
use Getopt::Long;
use Bio::SeqIO;
use File::Basename;
use Data::Dumper;

use Common;
use HMM;



#default paramters
#
my $trainFgFile = "";	#
my $trainBgFile = "";

#trainFgFile and trainBgFile was used to estimate the model
#so if the model file is given, the training data is not necessary
my $modelFile = "";		#pre-existing model file
my $modelFile2 = "";	#model file to save to

my $minSiteNum = 2;

#maximum distance to consider a clustering
my $maxDist = 50; #for fg
my $maxDist2 = 1000;  #for bg

my $binNum = 20; #number of bins to discretize conservation and accessibility score

my $exonTypeNum = 4; #intron/intergenic, coding, 5'utr, 3'utr
my $verbose = 0;
my $debug = 0;

my $noConservation = 0;
my $noRNAStructure = 0;
my $headToTail = 0;


##########################################################



GetOptions ('f:s'=>\$trainFgFile,
		'b:s'=>\$trainBgFile,
		'bin:i'=>\$binNum,
		'm|min-site:i'=>\$minSiteNum,
		'no-conservation'=>\$noConservation,
		'no-rna-structure'=>\$noRNAStructure,
		'head-to-tail'=>\$headToTail,
		'load-model:s'=>\$modelFile,
		'save-model:s'=>\$modelFile2,
		'd:i'=>\$maxDist,
		'v|verbose'=>\$verbose
);

my $prog = basename ($0);


my $getModel = 0;

if (((((-f $trainFgFile) && -f ($trainBgFile)) || (-f $modelFile))) && $modelFile2)
{
	$getModel = 1;		
}

if (@ARGV != 3 && !$getModel)
{
	print "identify conserved motif clusters using hidden markov model\n";
	print "Usage: $prog [options] <testSite.txt> <site.out.bed> <cluster.out.bed>\n";
	print "[OPTIONS]:\n";
	print " -f              [string]: foreground (postive) site data file\n";
	print " -b              [string]: background (negative) site data file\n";
	print " -bin            [int]   : bin number to discretize conservation/accessibility ($binNum)\n";
	print " --head-to-tail          : calculate the head-to-tail distance of neighbor sites (should not use if two sites can overlap)\n";
	print "                         : (default is head-to-head)\n";
	print " -m              [int]   : minimum number of sites in a cluster\n";
	print " -no-conservation        : do not use conservation info\n";
	print " -no-rna-structure       : do not use rna structure info\n";
	print " --load-model    [string]: model file to be loaded from\n";
	print " --save-model    [string]: modle file to be saved to\n";
	print " -d              [int]   : maximal distance allowed ($maxDist)\n";
	print " -v                      : verbose\n";
	exit (1);
}

$maxDist2 = $maxDist if $maxDist > $maxDist2;

my ($testFile, $outSiteFile, $outClusterFile) = @ARGV;


#//////////////////////////////////////////////////////////////////
#
# get Model parameters
#


my $model ={};
if (-f $modelFile)
{
	#read model
	print "reading model parameters from $modelFile ...\n" if $verbose;
	$model = readModel ($modelFile);	
	$binNum = $model->{"binNum"};
	$maxDist = $model->{"maxDist"};
	$maxDist2 = $model->{"maxDist2"};
	print "model loaded\n" if $verbose;
}
else
{
	Carp::croak "$trainFgFile does not exists\n" unless -f $trainFgFile;
	Carp::croak "$trainBgFile does not exists\n" unless -f $trainBgFile;

	$model->{"fg"} = estimateModelParam ($trainFgFile, $maxDist2, 0);
	$model->{"bg"} = estimateModelParam ($trainBgFile, $maxDist2, 0);

	for (my $i = $maxDist + 1; $i <= $maxDist2; $i++)
	{
		$model->{"fg"}->{"distance"}->[$i] = 0;
		$model->{"bg"}->{"distance"}->[$i] = 1e-16 if $model->{"bg"}->{"distance"}->[$i] <= 0;
	}

	$model->{"fg"}->{"mu"} = 3 if $model->{"fg"}->{"mu"} < 3;
	$model->{"bg"}->{"mu"} = $model->{"fg"}->{"mu"} * $model->{"bg"}->{"density"} / $model->{"fg"}->{"density"}; 

	if ($model->{"bg"}->{"mu"} < 2)
	{
		$model->{"bg"}->{"mu"} = 2;
		$model->{"fg"}->{"mu"} = $model->{"bg"}->{"mu"} * $model->{"fg"}->{"density"} / $model->{"bg"}->{"density"};
	}

	##NOTES:
	##this is not the best way, because it assumes that fg and bg seq have the same length
	
	##the following line gives an alternative
	#$model->{"bg"}->{"mu"} = $model->{"bg"}->{"density"} * 2000; #assume the average length of a bg seq is 5Kb, rough average length of introns
	Carp::croak "mu must be >=1: fg_mu=", $model->{"fg"}->{"mu"}, ", bg_mu=", $model->{"bg"}->{"mu"}, "\n" if $model->{"fg"}->{"mu"} < 1 || $model->{"bg"}->{"mu"} < 1;
}


if ($modelFile2)
{
	print "save HMM model to $modelFile2 ...\n" if $verbose;
	my $fout;
	
	open ($fout, ">$modelFile2") || Carp::croak "can not open file $modelFile2 to write\n";
	
	print $fout join ("\t", "bin_num", $binNum), "\n";
	print $fout join ("\t", "max_distance_positive", $maxDist), "\n";
	print $fout join ("\t", "max_distance_negative", $maxDist2), "\n";
	print $fout join ("\t", "mu", $model->{"fg"}->{"mu"}, $model->{"bg"}->{"mu"}), "\n";
	print $fout join ("\t", "distance_positive", @{$model->{"fg"}->{"distance"}}), "\n";
	print $fout join ("\t", "distance_negative", @{$model->{"bg"}->{"distance"}}), "\n";

	my @tmp = (0..$binNum);
	for (my $i = 0; $i < @tmp; $i++)
	{
		$tmp[$i] /= $binNum;
	}
	
	print $fout join ("\t", "#", @tmp), "\n";
	
	#positive conservation
	#0 --intron/intergenic 1--coding exon 2 --5'utr 3 -- 3'utr
	for (my $i = 0; $i < @{$model->{"fg"}->{"conservation"}}; $i++)
	{
		print $fout join("\t", "conservation_positive_$i", @{$model->{"fg"}->{"conservation"}->[$i]}), "\n";
	}
	
	#negative conservation
	for (my $i = 0; $i < @{$model->{"bg"}->{"conservation"}}; $i++)
	{
		print $fout join("\t", "conservation_negative_$i", @{$model->{"bg"}->{"conservation"}->[$i]}), "\n";
	}
	
	print $fout join("\t", "accessibility_positive", @{$model->{"fg"}->{"accessibility"}}), "\n";
	print $fout join("\t", "accessibility_negative", @{$model->{"bg"}->{"accessibility"}}), "\n";
	
	close ($fout);
	Carp::croak "mu must be >=1\n" if $model->{"fg"}->{"mu"} < 1 || $model->{"bg"}->{"mu"} < 1;
}

exit (0) unless $testFile && (-f $testFile);


#the state symbol of HMM
my @state = qw(S+ S- + - I+ I-);


#//////////////////////////////////////////////////////////////////
#
#now build transition and emission probability matrices
#

my $fgmu = $model->{"fg"}->{"mu"};
my $bgmu = $model->{"bg"}->{"mu"};

print "fg mu=$fgmu, bgmu = $bgmu\n" if $verbose;

#transition matrix
#---------------------------
my @A = (
		[0, 0, 1,           0,          0,              0],
		[0, 0, 0,           1-1/$bgmu,	1/$bgmu,       	0],
		[0, 0, 1- 1/$fgmu , 0,          1/$fgmu/$bgmu,	1/$fgmu * (1- 1/$bgmu)],
		[0, 0, 0,         	1-1/$bgmu,  1/$bgmu,        0],
		[0, 0, 1,           0,          0,              0],
		[0, 0, 0,          	1-1/$bgmu,  1/$bgmu,        0]);




#emission matrix of distance
#----------------------------
my @BDistanceStart; #distance distribution (dumbing) of start state

my $D = @{$model->{"fg"}->{"distance"}};

$BDistanceStart[0] = 1;

for (my $d = 1; $d < $D; $d++)
{
	$BDistanceStart[$d] = 0; #$bgDist->{"B"}->[$d]; #1e-30;
}

#distance distribution for each state
my @BDistance = (
		\@BDistanceStart, 
		\@BDistanceStart,
		$model->{"fg"}->{"distance"},
		$model->{"bg"}->{"distance"},
		$model->{"bg"}->{"distance"},
		$model->{"bg"}->{"distance"});


#emission matrix of conservation
#--------------------------------------
my @BConservation = (
		$model->{"fg"}->{"conservation"},
		$model->{"bg"}->{"conservation"},
		$model->{"fg"}->{"conservation"},
		$model->{"bg"}->{"conservation"},
		$model->{"fg"}->{"conservation"},
		$model->{"bg"}->{"conservation"});

#emission matrix of accessibility
#--------------------------------
my @BAccessibility = (
		$model->{"fg"}->{"accessibility"},
		$model->{"bg"}->{"accessibility"},
		$model->{"fg"}->{"accessibility"},
		$model->{"bg"}->{"accessibility"},
		$model->{"fg"}->{"accessibility"},
		$model->{"bg"}->{"accessibility"});
	

#initiation
#-----------------------------
my @pi = (0.5, 0.5, 0, 0, 0, 0);


my @alphabet;

my $N = @pi;
my $alphabetSize = $maxDist + 1;

if ($verbose)
{
	print "A=\n";
	for (my $i = 0; $i < $N; $i++)
	{
		print join ("\t", @{$A[$i]}), "\n";
	}

	print "BDistance=\n";
	for (my $i = 0; $i < $N; $i++)
	{
		print join ("\t", @{$BDistance[$i]}), "\n";
	}

	print "BConservation=\n";

	for (my $i = 0; $i < $N; $i++)
	{
		my $conservation = $BConservation[$i];
		for (my $exonic = 0; $exonic <@$conservation; $exonic++)
		{
			print join ("\t", @{$conservation->[$exonic]}), "\n";
		}
	}

	print "pi=", join ("\t", @pi), "\n";
}


print "reading test site from $testFile ...\n" if $verbose;


my $testSites = readSiteFile ($testFile);
	

print "posterior and viterbi decoding ...\n";

my $ret = HMM::posteriorDecode ($testSites, \@A, \&emission, \@pi);
my $ret2 = HMM::viterbi ($testSites, \@A, \&emission, \@pi);

for (my $s = 0; $s < @$testSites && $debug; $s++)
{

	my $seqId = $testSites->[$s]->[0]->{"seqId"};

	print "\nseqId = $seqId\n";

	print "dist:";
	for (my $i = 0; $i < @{$testSites->[$s]}; $i++)
	{
		print "\t", $testSites->[$s]->[$i]->{"distance"};
	}
	print "\n";
	
	print "start:";
	for (my $i = 0; $i < @{$testSites->[$s]}; $i++)
	{
		print "\t", $testSites->[$s]->[$i]->{"start"};
	}

	print "\n";
	print "cons:";
	for (my $i = 0; $i < @{$testSites->[$s]}; $i++)
	{
		printf ("\t%.2f", $testSites->[$s]->[$i]->{"conservation"});
	}
	print "\n";
	
	print "exonic:";
	for (my $i = 0; $i < @{$testSites->[$s]}; $i++)
	{
		print "\t", $testSites->[$s]->[$i]->{"exonic"};
	}
	print "\n";
	
	print "access:";
	for (my $i = 0; $i < @{$testSites->[$s]}; $i++)
	{
		printf ("\t%.2f", $testSites->[$s]->[$i]->{"accessibility"});
	}
	print "\n";
	

	#print join ("\t", @{$testObs[$s]}), "\n";
	print "Prob:";
	
	for (my $i = 0; $i < @{$ret->[$s][0]}; $i++)
	{
		#my $state = $ret2->[$s]->[$i];
		my $score = $ret->[$s][0][$i] + $ret->[$s][2][$i] + $ret->[$s][4][$i];
		#printf ("\t%.1f", $ret->[$s][$state][$i]);
		printf ("\t%.2f", $score);
		#$ret->[$s][2][$i] = sprintf ("%.1f", $ret->[$s][2][$i]);
	}
	print "\n";
	#print "Prob:\t", join ("\t", @{$ret->[$s][2]}), "\n";
	print "State:\t", join ("\t", @state[@{$ret2->[$s]}]), "\n";
}

=test
print "dump individual sites into $outSiteFile ...\n";

my $fout = new FileHandle;
open ($fout, ">$outSiteFile") || Carp::croak "can not open file $outSiteFile to write\n";

print $fout "track name=sites useScore = 1\n";
for (my $s = 0; $s < @$testSites; $s++)
{
	my $sitesInSeq = $testSites->[$s];

	my $seqId = $sitesInSeq->[0]->{"seqId"};

	for (my $i = 0; $i < @$sitesInSeq; $i++)
	{
		my $site = $sitesInSeq->[$i];

		my $stateIdx = $ret2->[$s][$i];
		next unless $stateIdx == 0 || $stateIdx == 2 || $stateIdx == 4;
		
		my $score = ($ret->[$s][0][$i] + $ret->[$s][2][$i] + $ret->[$s][4][$i]) * 1000;
		print $fout join ("\t", $seqId, $site->{"start"}, $site->{"end"} + 1, $site->{"id"}, $score, "+"), "\n";
	}
}
close ($fout);
=cut



print "dump clusters and sites ...\n";

my $fout;
my $foutSite;
open ($fout, ">$outClusterFile") || Carp::croak "can not open file $outClusterFile to write\n";
print $fout "track name=cluster\n";

open ($foutSite, ">$outSiteFile") || Carp::croak "can not open file $outSiteFile to write\n";
print $foutSite "track name=sites useScore = 1\n";


for (my $s= 0; $s < @$testSites; $s++)
{
	my $sitesInSeq = $testSites->[$s];
	my $seqId = $sitesInSeq->[0]->{"seqId"};

	my @cluster;
	my $clusterScore = 0;
	my $clusterIter = 0;

	for (my $i = 0; $i < @$sitesInSeq; $i++)
	{
		my $stateIdx = $ret2->[$s][$i];
		my $site = $sitesInSeq->[$i];
		my $score = $ret->[$s][0][$i] + $ret->[$s][2][$i] + $ret->[$s][4][$i];
		
		if ($stateIdx != 0 && $stateIdx != 2) # not S+ or +, end a cluster, start a new cluster
		{
			if (@cluster >= $minSiteNum)
			{
				my $chromStart = $cluster[0]->{"start"};
				my $chromEnd = $cluster[$#cluster]->{"end"} + 1;
				#$clusterScore /= @cluster;
				my $clusterId = $seqId ."_$clusterIter" . sprintf ("[%.2f]", $clusterScore);

				print $fout join ("\t", $seqId, $chromStart, $chromEnd, $clusterId, $clusterScore, '+'), "\n";

				foreach my $siteInCluster (@cluster)
				{
					print $foutSite join ("\t", $seqId, $siteInCluster->{"start"}, $siteInCluster->{"end"} + 1, 
							$siteInCluster->{"id"} . "@" . $clusterId , $siteInCluster->{'score'} * 1000, "+"), "\n";
				}
				$clusterIter++;
			}
			@cluster = ();
			$clusterScore = 0;
		}
		
		if ($stateIdx == 0 || $stateIdx == 2 || $stateIdx == 4) #S+, + or I+, extend a cluster
		{
			my $score = $ret->[$s][0][$i] + $ret->[$s][2][$i] + $ret->[$s][4][$i];
			$site->{'score'} = $score;
			push @cluster, $site;
			#$clusterScore += log($score/(1-$score+0.0000000000001)) / log(2);

			$clusterScore += log(emission ($stateIdx, $site) / emission ($stateIdx+1, $site));
		}
	}

	if (@cluster >= $minSiteNum)
	{
		my $chromStart = $cluster[0]->{"start"};
		my $chromEnd = $cluster[$#cluster]->{"end"} + 1;
		#$clusterScore /= @cluster;
		my $clusterId = $seqId ."_$clusterIter" . sprintf ("[%.2f]", $clusterScore);

		print $fout join ("\t", $seqId, $chromStart, $chromEnd, $clusterId, $clusterScore, '+'), "\n";

		foreach my $siteInCluster (@cluster)
		{
			print $foutSite join ("\t", $seqId, $siteInCluster->{"start"}, $siteInCluster->{"end"} + 1, 
					$siteInCluster->{"id"} . "@" . $clusterId , $siteInCluster->{'score'} * 1000, "+"), "\n";
		}
	}
}

close ($fout);
close ($foutSite);


sub readModel
{
	my $modelFile = $_[0];
	my $fin;

	open ($fin, "<$modelFile") || Carp::croak "cannot open file $modelFile to read\n";
	
	my $model = {};
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/\#/;
		my @values = split (/\t/, $line);
		my $name = shift @values;
		if ($name eq "mu")
		{
			$model->{"fg"}->{"mu"} = $values[0];
			$model->{"bg"}->{"mu"} = $values[1];
		}		
		elsif ($name eq "bin_num")
		{
			$model->{"binNum"} = $values[0];
		}
		elsif ($name eq 'max_distance_positive')
		{
			$model->{"maxDist"} = $values[0];
		}
		elsif ($name eq 'max_distance_negative')
		{
			$model->{"maxDist2"} = $values[0];
		}
		elsif ($name eq 'distance_positive')
		{
			$model->{"fg"}->{"distance"} = \@values;
		}
		elsif ($name eq 'distance_negative')
		{
			$model->{"bg"}->{"distance"} = \@values;
		}
		elsif ($name =~/^conservation_positive_(\d+)$/)
		{
			my $i = $1;
			$model->{"fg"}->{"conservation"}->[$i] = \@values;
		}
		elsif ($name =~/^conservation_negative_(\d+)$/)
		{
			my $i = $1;
			$model->{"bg"}->{"conservation"}->[$i] = \@values;
		}
		elsif ($name eq 'accessibility_positive')
		{
			$model->{"fg"}->{"accessibility"} = \@values;
		}
		elsif ($name eq 'accessibility_negative')
		{
			$model->{"bg"}->{"accessibility"} = \@values;
		}
		else
		{
			Carp::croak "incorrect parameter : $name\n";
		}
	}
	close ($fin);
	return $model;
}

sub emission
{
	my ($state, $obs) = @_;
	my $d = $obs->{"distance"};
	#return 1e-15 if $d > $maxDist;
	$d = $maxDist if $d > $maxDist2;
	my $c = $obs->{"conservation"};
	my $a = $obs->{"accessibility"};
	my $exonic = $obs->{"exonic"};

	my $idxC = discretize ($c, 0, 1, $binNum);
	my $idxA = discretize ($a, 0, 1, $binNum);

	#print "state = $state, d = $d, c = $c\n";
	#exonic: 0 -- intronic/intergenic, 1-coding exon, 2- 5'utr exon, 3- 3'utr exon
	if ($noConservation && $noRNAStructure)
	{
		return $BDistance[$state][$d];
	}
	elsif ($noConservation)
	{
		return $BDistance[$state][$d] * $BAccessibility[$state][$idxA];
	}
	elsif ($noRNAStructure)
	{
		return $BDistance[$state][$d] * $BConservation[$state][$exonic][$idxC];
	}
	else
	{
		return $BDistance[$state][$d] * $BConservation[$state][$exonic][$idxC] * $BAccessibility[$state][$idxA];
	}
}

sub discretize
{
	my ($v, $min, $max, $n) = @_;
	
	my $binSize = ($max - $min) / $n;
	return int(($v-$min)/$binSize+0.5);
}

#if truncate, find truncated distribution
#otherwise, collpase all distances beyond the limit
sub estimateModelParam
{
	my ($siteFile, $maxDistance, $truncate) = @_;

	print "reading sites from $siteFile ...\n" if $verbose;
	my $sites = readSiteFile ($siteFile);
	
	my @distance;
	my @conservation;

	#my @conservationExonic;
	#my @conservationIntronic;
	my @accessibility;
	my $totalSeqLen = 0;

	for (my $i = 0; $i <= $maxDistance; $i++)
	{
		$distance[$i] = 1e-15;
	}
	
	for (my $i = 0; $i <= $binNum; $i++)
	{
		for (my $exonic = 0; $exonic < $exonTypeNum; $exonic++)
		{
			$conservation[$exonic][$i] = 1e-15;
		}
		$accessibility[$i] = 1e-15;
	}
	
	foreach my $sitesInSeq (@$sites)
	{

		my $lastIdx = @$sitesInSeq - 1;
		$totalSeqLen += $sitesInSeq->[$lastIdx]->{"end"} - $sitesInSeq->[0]->{"start"} + 1;

		foreach my $site (@$sitesInSeq)
		{
			my $s = $site->{"start"};
			my $c = $site->{"conservation"};
			my $a = $site->{"accessibility"};
			my $exonic = $site->{"exonic"};

			my $d = $site->{"distance"};
			
			if ($d <= $maxDistance)
			{
				$distance[$d]+=1 if $d > 0;
			}
			else
			{
				$distance[$maxDistance]+=1 if $truncate == 0 && $d > 0;
			}

			my $idx = discretize ($c, 0, 1, $binNum);
		
			$conservation[$exonic][$idx]+=1;	
			
			$idx = discretize ($a, 0, 1, $binNum);
			$accessibility[$idx]+=1;
		}
	}

	#scaling
	my $sum = sum (\@distance);

	Carp::croak "one site in each sequence only?\n" if $sum < 2;

	print int($sum),  " pairs of sites with distance <= $maxDistance\n" if $verbose;

	for (my $i = 0; $i < @distance; $i++)
	{
		$distance[$i] /= $sum;
	}

	for (my $exonic = 0; $exonic < $exonTypeNum; $exonic++)
	{
		$sum = sum ($conservation[$exonic]);
	
		print int($sum), " sites of type $exonic\n" if $verbose;

		for (my $i = 0; $i < @{$conservation[$exonic]}; $i++)
		{
			$conservation[$exonic][$i] /= $sum;
		}
	}

	$sum = sum (\@accessibility);

	print int($sum), " sites with accessibility scores\n" if $verbose;

	for (my $i = 0; $i < @accessibility; $i++)
	{
		$accessibility[$i] /= $sum;
	}
	
	my $seqNum = @$sites;	
	
	my $avgSiteNumPerSeq = $sum / $seqNum;
	my $avgSiteDensity = $sum / $totalSeqLen;

	return {distance=>\@distance, conservation=>\@conservation, 
			accessibility=>\@accessibility, mu=>$avgSiteNumPerSeq, density=>$avgSiteDensity};
}

sub readSiteFile
{
	my $inFile = $_[0];
	my %sites;
	my @sites;

	my $fin;
	open ($fin, "<$inFile") || Carp::croak "cannot open file $inFile\n";
	my $i = 0;

	my $prevSeqId = "";
	my $prevSiteStart = -1;

	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;

		print "$i...\n" if $verbose && $i % 50000 == 0;
		$i++;

		my ($seqId, $siteStart, $siteEnd, $siteId, $score, $strand, $region, $bls, $pu, $exonic) = split (/\t/, $line);
		$exonic = 2 if $exonic == 4; #treat lnRNA exon as 3'UTR

		$siteEnd--;

		my $siteDist = 0;
		push @{$sites{$seqId}}, {seqId=> $seqId, id=>$siteId, start=>$siteStart, end=>$siteEnd, conservation=>$bls, accessibility=>$pu, exonic=>$exonic};
	}
	close ($fin);

	#sort sites and calculate the distance between sites; the first site in each seq have a dummy number of zero
	foreach my $seqId (keys %sites)
	{
		my $sitesInSeq = $sites{$seqId};

		my @sitesInSeqSorted = sort {$a->{"end"} <=> $b->{"end"}} @$sitesInSeq;
		@sitesInSeqSorted = sort {$a->{"start"} <=> $b->{"start"}} @sitesInSeqSorted;
		
		my $sPrev = -1;

		$sitesInSeqSorted[0]->{"distance"} = 0;

		for (my $i = 1; $i < @sitesInSeqSorted; $i++)
		{
			$sitesInSeqSorted[$i]->{"distance"} = $headToTail ? ($sitesInSeqSorted[$i]->{"start"} - $sitesInSeqSorted[$i-1]->{"end"}) 
				: ($sitesInSeqSorted[$i]->{"start"} - $sitesInSeqSorted[$i-1]->{"start"});
			Carp::croak "negative distance detected\n" if $sitesInSeqSorted[$i]->{"distance"} <= 0;
		}
		push @sites, \@sitesInSeqSorted;
	}
	return (\@sites);
}


