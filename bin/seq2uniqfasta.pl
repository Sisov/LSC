#!/usr/bin/perl -w
# Pre: Take as STDIN the sequences, and as arugments:
#      the outputfile
#      the temporary directory
#      the script diretory
# Post: Write the unique fasta file
use strict;
if(scalar(@ARGV) != 3) { die; }
my $outfile = shift @ARGV;
my $tdir = shift @ARGV;
my $scriptdir = shift @ARGV;
open(OF,"| sort -T $tdir | uniq -c | ".$scriptdir."uniqcount2fasta.pl > $outfile") or die;
while(my $line = <STDIN>) {
  print OF $line;
}
close OF;
