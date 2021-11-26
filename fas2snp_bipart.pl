#!/usr/bin/env perl -w
#

#die "usage: fas2snp_bipart.pl <FASTA> \n" if (!defined(@ARGV));
die "usage: fas2snp_bipart.pl <FASTA> \n" if ( $#ARGV != 0 );

chomp($fasta_file = $ARGV[0]);

open(FILE1, "<$fasta_file" ) or die "Can't open $fasta_file : $!";
chomp(@FASTA = <FILE1>);
close(FILE1);


%aGen = ();

foreach $line (@FASTA) {
    if ($line =~ />(\S+)/) {
	$currSeq = $1;
	$aGen{$currSeq} = "";
    }else {
	$line =~s/\s//g;
	$line =~tr/acgtn/ACGTN/;
	$aGen{$currSeq} .= $line;
    }
}

#create work hash with snp data 
    
foreach $sq ( sort {$a cmp $b } keys(%aGen)) {
    
#    $aGen{$sq} =~ tr/atcgn/ATCGN/;
#    $aGen{$sq} =~ s/\s//g;

    @tmpArray = split(//, $aGen{$sq});
    for ($i = 0 ; $i <= $#tmpArray; $i++) {
	$tmpHash{$i}{$sq}= $tmpArray[$i];
    }
}
   

foreach $seqpos ( sort { $a <=> $b } keys(%tmpHash)) {

    %tmpBases1 = ();
    %tmpBases2 = ();
    %tmpSeqs = ();

    foreach $sq ( sort { $a cmp $b } keys(%{$tmpHash{$seqpos}})) {
	if($tmpHash{$seqpos}{$sq} ne "N" && $tmpHash{$seqpos}{$sq} ne "-") {
	    $tmpBases1{$tmpHash{$seqpos}{$sq}}{$sq} = 1;
	}
	$tmpBases2{$tmpHash{$seqpos}{$sq}} = 1;
	$tmpSeqs{$sq}=$tmpHash{$seqpos}{$sq};
    }
    $noVariants = scalar keys %tmpBases1;
    if ($noVariants == 2) {
	@bipart = ();
	@sbipart = ();
	foreach $base ( sort { $a cmp $b } keys %tmpBases1 ) {
	    @tmpArray = ();
	    @stmpArray = ();
	    $tmpstr = "";
	    
	    foreach $sq ( sort { $a cmp $b } keys %{$tmpBases1{$base}} ) {
		push @tmpArray, $sq;
	    }
	    @stmpArray = sort @tmpArray;
	    
	    $tmpstr =  join(':', @stmpArray);
	    
	    push @bipart, $tmpstr;
	}
	@sbipart = sort @bipart;
	
	@s1 = split(/:/, $sbipart[0]);
	@s2 = split(/:/, $sbipart[1]);

	foreach $strain (sort { $a cmp $b } keys(%tmpSeqs)){
	    print $tmpSeqs{$strain} . "\t";
	}
	print "$seqpos\t$sbipart[0]\t$tmpSeqs{$s1[0]}\t$sbipart[1]\t$tmpSeqs{$s2[0]}\n";
	
    }
}



