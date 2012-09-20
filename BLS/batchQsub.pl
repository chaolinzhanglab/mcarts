#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Carp;
use File::Basename;

use MyConfig;
use SGE;


my $prog = basename ($0);


my $cache = MyConfig::getDefaultCache ($prog); 

my $jobName = $ENV{'LOGNAME'};
my $queueFile = "";
my $dir = "";
my $qsubOK = 0;
my $memory = "";
my $noSGE = 0;
my $waitForSGE = 0;
my $verbose = 0;

GetOptions ('c|cache:s'=>\$cache,
		'j|jobname:s'=>\$jobName,
		'm|memory:s'=>\$memory,
		'q|queues:s'=>\$queueFile,
		'd|dir:s'=>\$dir,
		'wait'=>\$waitForSGE,
		'no-sge'=>\$noSGE,
		'v'=>\$verbose);


if (@ARGV != 1)
{
	print "submit jobs via qsub\n";
	print "Usage: $prog [options] <script.list>\n";
	print " -c [string]: cache dir for qsub output [$cache]\n";
	print " -j [string]: job names [$jobName]\n";
	print " -m [string]: memory (e.g. 2G, default is not to specify)\n";
	print " -d [string]: dir to be attached to each script [$dir]\n";
	print " -q [string]: file name specifying queques to be used\n";
	print " --wait     : wait until the jobs are finished before exiting\n";
	print " --no-sge   : do not submit jobs to SGE; run locally\n";
	print " -v         : verbose\n";
	exit(1);
}


die "$cache already exist\n" if (-d $cache);
die "$cache already exist\n" if (-f $cache);

my $queueStr = "";

if ($noSGE == 0)
{
	my $testQsub = `which qsub`;
	chomp $testQsub;
	$qsubOK = 1 if ($testQsub =~/\/qsub$/);
}

system ("mkdir $cache") if $qsubOK;

if (-f $queueFile)
{
	print "reading queue information from $queueFile ...\n" if $verbose;
	my $fin;
	open ($fin, "<$queueFile") || Carp::croak "can not open file $queueFile to read\n";
	my @queues =<$fin>;
	close ($fin);

	chomp @queues;
	$queueStr = join (",", @queues);
}

my $scriptListFile = $ARGV [0];

my $fin;

my @jobIds;

open ($fin, "<$scriptListFile") || Carp::croak "can not open file $scriptListFile to read\n";
while (my $f = <$fin>)
{
	chomp $f;
	next if $f=~/^\s*$/;
	
	if (length($dir) > 0 && (-d $dir))
	{
		$f = "$dir/$f";
	}
	my $base = basename ($f);
	
	my $memOpt = "";
	$memOpt = "-l virtual_free=$memory" if $memory;

	#print "qsub -cwd -e $cache -o $cache -N $jobName $f\n";
	if ($qsubOK)
	{

		#use default queues
		my $queueFlag = $queueStr ne '' ? '-q $queueStr' : '';

		#system ("qsub -cwd $memOpt -v TMPDIR=$cache -e $cache -o $cache -N $jobName.$base $f");
		#system ("qsub -cwd $memOpt -v TMPDIR=$cache -V -e $cache -o $cache -N $jobName.$base $f");
		my $ret = `qsub -cwd $memOpt -v $queueFlag TMPDIR=$cache -V -e $cache -o $cache -N $jobName.$base $f`;
		chomp $ret;
		my @cols = split (/\s+/, $ret);
		my $jobId = $cols[2];
		Carp::croak "invalid job id: $jobId\n" unless $jobId =~/\d+/;

		push @jobIds, $jobId;
		
		print "job id=$jobId submited ...\n" if $verbose;
		#system ("/chome/sge/bin/glinux/qsub -cwd -v TMPDIR=$cache -V -e $cache -o $cache -N $jobName $f");
	}
	else
	{
		print "sh $f> $f.out 2> $f.err\n";
		my $ret = system ("sh $f > /dev/null");
		
		#my $ret = system ("sh $f > $f.out 2> $f.err");
		Carp::croak "jobs failed: $!\n" if $ret != 0;
	}
}
close ($fin);

if ($qsubOK && $waitForSGE && @jobIds > 0)
{
	print "job ids=", join ("\t", @jobIds), "\n" if $verbose;
	waitUntilSGEJobsDone (\@jobIds, $verbose);
}

