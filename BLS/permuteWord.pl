#!/usr/bin/perl -w

use strict;
use Common;


die "usage: script <consensus>\n" if @ARGV < 1;
my $motif = $ARGV[0];
my $len = length ($motif);
my $randIdx = randSeq (0, $len);
my @tmp = split (//, $motif);
my @tmp2 = @tmp [@$randIdx];
my $motif_permute = join ("", @tmp2);
print $motif_permute, "\n";


