#!/usr/bin/perl -w


use strict;
use warnings;

use File::Basename;

use Carp;
use Data::Dumper;
use Getopt::Long;

use PhyloTree;


my $prog = basename ($0);
my $debug = 0;
my $verbose = 0;

GetOptions (
		"v|verbose"=>\$verbose);

if (@ARGV != 3)
{
	print "evaluate branch length score of motif sites\n";
	print "Usage: $prog [options] <tree.nh> <site.tab> <site.out>\n";
	print "OPTIONS\n";
	print " -v : verbose\n";
	exit (1);
}


my ($treeFile, $siteFile, $siteConsFile) = @ARGV;


#-----------------------------------------------------------------------
#load phylogenetic tree

print "load phylogentic tree from $treeFile ...\n" if $verbose;

my $tree = readTreeFile ($treeFile);

Carp::croak "empty tree?\n" if $tree == 0;

my $leaves = getLeafNodes ($tree);

my $nspec = keys %$leaves;

print "the tree has $nspec species\n" if $verbose;

Carp::croak "too few species\n" if $nspec < 2;



#----------------------------------------------------------------------
#load site
#
print "loading site information from $siteFile ...\n" if $verbose;

my @sites;

my $fin;
open ($fin, "<$siteFile") || Carp::croak "can not open file $siteFile\n";

#read the species information
my $header = <$fin>;
chomp $header;
my @species = split ("\t", $header);
shift @species;

my %speciesHash;
foreach my $s (@species)
{
	Carp::croak "species $s is not included in the tree\n" unless exists $leaves->{$s};
	$speciesHash{$s} = 1;
}

while (my $line =<$fin>)
{
	chomp $line;
	my @cols = split ("\t", $line);
	push @sites, \@cols;
}

close ($fin);

my $nsites = @sites;

print "$nsites sites loaded\n" if $verbose;


my @speciesToRemove;
#prune irrelevant species in the tree
foreach my $s (keys %$leaves)
{
	push @speciesToRemove, $leaves->{$s} unless exists $speciesHash{$s};
}

$tree = removeLeafNodes ($tree, \@speciesToRemove);
#$leaves = Common::getLeafNodes ($tree);
my $totalBL = totalBranchLength ($tree);


print "the total branch length of the tree is $totalBL\n" if $verbose;


print "calculating branch length score for each site\n" if $verbose;
my $fout;

open ($fout, ">$siteConsFile") || Carp::croak "can not open file $siteConsFile ...\n";

my $iter = 0;
print $fout join ("\t", "Site", @species, "BLS"), "\n";
foreach my $s (@sites)
{
	$iter++;

	if ($iter - int($iter/1000) * 1000 == 0)
	{
		print $iter, "...\n" if $verbose;
	}
	
	my @oneSite = @$s;
	my $name = shift @oneSite;
	my @speciesWithoutSite;

	my $ncol = @oneSite;
	#print join("\t", $name, @oneSite), "\n";
	#print "sum=", Common::sum (\@oneSite), ", ncol = $ncol\n";
	
	for (my $i = 0; $i < @oneSite; $i++)
	{
		push @speciesWithoutSite, $species[$i] if $oneSite[$i] == 0;
	}
	#my $n = @speciesWithoutSite;
	#print "$n species without site:\n" if $verbose;

	my $treeCpy = copyTree ($tree);
	my $siteTree = removeSpecies ($treeCpy, \@speciesWithoutSite);
	
	my $siteBL = totalBranchLength ($siteTree);
	
	releaseTree ($siteTree);
	
	#print "tree=", Common::codeTree ($siteTree), ", blen = $siteBL\n";
	
	my $BLS = $siteBL / $totalBL;

	print $fout join ("\t", $name, @oneSite, $BLS), "\n";
}

close ($fout);




