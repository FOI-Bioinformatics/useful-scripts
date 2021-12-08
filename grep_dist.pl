#!/usr/bin/env/ perl -w
#
#	split a multifasta file to separate fasta files 
#
#
use List::Util qw/ min max sum /;

die "usage: grep_dist.pl <distances> <refs> <dist>\n" if (!defined(@ARGV));
die "usage: grep_dist.pl <distances> <refs> <dist>\n" if ( $#ARGV != 2 );

chomp($dst = $ARGV[0]);
chomp($rf = $ARGV[1]);
chomp($dist = $ARGV[2]);

open(FILE1, "<$dst" ) or die "Can't open $dst : $!";
chomp(@A = <FILE1>);
close(FILE1);

open(FILE1, "<$rf" ) or die "Can't open $rf : $!";
chomp(@RF = <FILE1>);
close(FILE1);

foreach $line (@A) {
    
    @tmp = split(/\s+/, $line);
    $dists{$tmp[0]}{$tmp[1]}=$tmp[2];
    $dists{$tmp[1]}{$tmp[0]}=$tmp[2];    
}

foreach $st1 (@RF) {
    foreach $st2 ( sort {$a cmp $b } keys(%{$dists{$st1}})) {
	if($dists{$st1}{$st2} <= $dist) {
	    $out{$st1} = 1;
	    $out{$st2} = 1;
	}
    }
}

foreach $sq ( sort {$a cmp $b } keys(%out)) {
    print "$sq \n";
}
