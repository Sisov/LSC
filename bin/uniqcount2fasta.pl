#!/usr/bin/perl -w
use strict;

# Pre read the input of sequences from stdin of a uniq -c of sequences
# Post output a fasta file of unique sequences

my $z = 0;
while(my $line = <STDIN>) {
  chomp($line);
  $z += 1;
  if($line=~/^\s*(\d+)\s+(\S+)$/) {
    print ">$z".'_'.$1."\n".$2."\n";
  }
}
