#!/usr/bin/perl -w

=head1 Name

    SISP.pl
    Version: 0.1
    Date: June.7 2013
    Contact: kevinchjp@gmail.com

=head1 Description

    This program is used to quickly identify strain species identity

=head1 Usage

    perl SISP.pl -bl blastn -fd formatdb -gn genome -od output directory
    Arguments explained
    bl: Directory of blastn excecutable file
    fd: Directory of BLAST formatdb excecutable file
    gn: genome sequence file in FASTA format for query strain
    od: output directory

=cut

use strict;
use FindBin qw($Bin);
use Getopt::Long;

my ($bl,$fd,$gn,$od);
GetOptions(
    "bl=s" => \$bl,
    "fd=s" => \$fd,
    "gn=s" => \$gn,
    "od=s" => \$od
          );

die `pod2text $0` unless $bl && $fd && $gn && $od;

#split query genome to speed up BLAST alignment
my $queryGsize; my $chop_len = 1020;
unless(-d $od){`mkdir $od`;}
$/ = "\n>";
open SR,$gn or die "$gn $!\n";
open CR,">$od/ref.split";
open CL,">$od/ref.split.len";
while(<SR>){
    chomp;
    s/>//g;
    my ($scaf,$seq) = split /\n/,$_,2;
    my $scaf_name = (split /\s+/,$scaf)[0];
    $seq =~ s/\s+//g;
    my $seq_len = length($seq);
    $queryGsize += $seq_len;
    my @cut=($seq =~ /(.{1,$chop_len})/g);
    for my $cut_num (0..$#cut){
        next unless $cut[$cut_num] =~ /[ATCGatcg]+/;
        my $cut_len = length($cut[$cut_num]);
        next if $cut_len < 100; 
        print CR ">$scaf_name\_$cut_num\n$cut[$cut_num]\n";
        print CL "$scaf_name\_$cut_num\t$cut_len\n";
    }
}
close SR;close CR;close CL;
$/ = "\n";

#read Hierarchy.db and store information in memory
open IH,"$Bin/Hierarchy.db" or die "Hierarchy.db $!\n";
my %system;
while(<IH>){
    chomp; my @t=split;
    push @{$system{$t[0]}{$t[1]}{$t[2]}} , $t[3];
}
close IH;

#read RS.tab and store information in memory
open IN,"<$Bin/RS.tab" or die "RS.tab $!\n";
my %rs;
while(<IN>){
    chomp;
    my @t=split /\s+/, $_;
    $rs{$t[0]} = $t[1];
}
close IN;

my $count=1; my %str_tni; our $maxTni=0; our $maxStr="initial"; #$count: how many db strains have been compared; %str_tni: $str_tni{Strain genome dir}=Its TNI with query; $maxTni: current maximum TNI; $maxStr: strain with max TNI
foreach my $flg (keys %system){ #$flg:first level groups ID
    my $tni = &TNI($count,$rs{$flg});
    &Max($tni,$rs{$flg});
    #go in the group
    if ($tni > 3.2){
        foreach my $slg (keys %{$system{$flg}}){
            $count++;
            my $tni = &TNI($count,$rs{$slg});
            &Max($tni,$rs{$slg});
            #go to the second level
            if ($tni > 12.8){
                foreach my $tlg (keys %{$system{$flg}{$slg}}){
                    $count++;
                    my $tnirs = &TNI($count,$rs{$tlg});
                    &Max($tnirs,$rs{$tlg});
                    next if $tnirs < 51.2;
                    #go to the third level
                    foreach my $tlgs (@{$system{$flg}{$slg}{$tlg}}) { #$tlgs: third level group strains
                        next if $tlgs eq $tnirs; #next if current str is representative str
                        $count++;
                        my $tni = &TNI($count,$tlgs);
                        &Max($tni,$tlgs);
                    }
                    last;
                }
                last;
            }
        }
        last;
    }
}

print "Maximum TNI strain: $maxStr\nTNI: $maxTni\n";
#subrouting to update max TNI and strain report
sub Max{
    my ($tni,$str) = @_;
    $str_tni{$str}=$tni;
    if ($tni>$maxTni){$maxTni=$tni; $maxStr=$str;}
    print "$str\t$tni\n";
}

#subrouting for calculating TNI
sub TNI{
    my ($n,$pathway) = @_;
    `ln -sf $Bin/$pathway $od/s$n.fa`;
    `$fd -i $od/s$n.fa -p F`;
    $/ = "\n>";
    my ($tagGsize,$small_size); #tagGsize: genome size for target strain (db strain); small_size: smaller size of query and target
    open SR,"$od/s$n.fa" or die "s$n.fa $!\n";
    open SL,">$od/s$n.len";
    while(<SR>){
        chomp;
        s/>//g;
        my ($scaf,$seq) = split /\n/,$_,2;
        my $scaf_name = (split /\s+/,$scaf)[0];
        $seq =~ s/\s+//g;
        my $seq_len = length($seq);
        $tagGsize += $seq_len;
        print SL "$scaf_name\t$seq_len\n";
    }
    $/ = "\n";
    close SR;close SL;
    if($queryGsize > $tagGsize){
        $small_size = $tagGsize;
    }
    else{$small_size = $queryGsize;}
    `$bl -i $od/ref.split -d $od/s$n.fa -X 150 -q -1 -F F -e 1e-15 -m 8 -a 2 -o $od/s$n.blast -p blastn`;
    my $index = `$Bin/Blast_filter.pl $small_size $od/s$n.blast $od/s$n.len $od/ref.split.len`;
    `rm $od/s$n.fa* $od/s$n.blast $od/s$n.len`;
    return $index;
}
