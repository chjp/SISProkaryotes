#! /usr/bin/perl -w
use strict;

die "This program is to filter the blast and get the taxonomy index.
Usage:  perl $0 <smaller genome size> <blastout> <dbase length> <query length>\n" unless @ARGV ==4;

my $identity = 30;
my $cvg = 70;
my (%dup,%dlen,%qlen);
open DL,$ARGV[2] or die "$ARGV[2],Open Wrong!\n";
while(<DL>){
    chomp;
    my @t = split;
    $dlen{$t[0]} = $t[1];
}
close DL;
open QL,$ARGV[3] or die "$ARGV[3],Open Wrong!\n";
while(<QL>){
    chomp;
    my @t = split;
    $qlen{$t[0]} = $t[1];
}
close QL;

open BT,$ARGV[1] or die "$ARGV[1],Open Wrong!\n";
my $sum = 0;my $count = 0;
my %subdup;
my $subjectdup = 0;
my $mapLen = 0;
my $matchNum = 0;
my $Hr = 0;
my $size = $ARGV[0];
my ($ANI,$HrRatio,$mapLenRatio,$matchNumRatio);
while(<BT>){
    chomp;
    my @t = split;
    next if $t[3] < 100;
    unless($#t == 11){
        die "$ARGV[1] uncomplete!\n";
    }
    unless( (exists $dlen{$t[1]}) && (exists $qlen{$t[0]}) ){
        die "$ARGV[1] wrong format!\n";
    }
    next if exists $dup{$t[0]};
    $dup{$t[0]} = 1;
    next if $t[2]<=$identity;
    my $t1 = $t[3]*100/$dlen{$t[1]};
    my $t2 = $t[3]*100/$qlen{$t[0]};
    if($t2 > $t1){$t1 = $t2;}
    next if $t1<$cvg;
    # $sum += $t[2];
    #$count++;

    my $subflag = 0;
    if($t[8] > $t[9]){
        my $temp = $t[8];
        $t[8] = $t[9];
        $t[9] = $temp;
    }
    if(exists $subdup{$t[1]}){
        for my $num (keys %{$subdup{$t[1]}} ){
            my ($x1,$x2) = split /\t/,$subdup{$t[1]}{$num};
            if( (($t[8] >= $x1) && ($t[8] <= $x2)) || (($t[9] >= $x1) && ($t[9] <= $x2)) ){
                $subflag = 1;
            }
        }
    }
    next if $subflag == 1;
    $subdup{$t[1]}{$subjectdup++} = "$t[8]\t$t[9]";

    $sum += $t[2];
    $count++;
    $Hr += $qlen{$t[0]};
    $mapLen += $t[3];
    $matchNum += $t[3]*$t[2];
}
close BT;
$HrRatio = $Hr*100/$size;
$mapLenRatio = $mapLen*100/$size;
$matchNumRatio = $matchNum/$size;
if($mapLenRatio >100){
    $mapLenRatio = 100;
}
if($matchNumRatio>100){
    $matchNumRatio = 100;
}
if($count == 0){
    $ANI = 30;
}
else{$ANI = $sum/$count;}
print "$matchNumRatio";
