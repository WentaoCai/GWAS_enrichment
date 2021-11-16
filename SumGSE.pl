#!/usr/bin/perl
###########################################################################
## SumGSE --- Process genome region bed, GWAS summary
############################################################################
## Author: Wentao Cai <wtaocai@gmail.com>
#use warnings;
#use strict;
print " 

#####################################
#                                   #
# SumGSE v0.02                      #
#                                   #
# last change: 17/10/2021           #
# By Wentao Cai <wtaocai\@gmail.com> #
#####################################

";

my $usage = "\nUsage:\n\t perl SumGSE.pl -i [genome_region.bed] -g [GWAS_summaries.txt] -s [specific.regions] -e [genome_region_extention] -n [permutation_times]\n;


Options:

        -i    input file in bed format. The input file should be genome feature regions(such as DEGs，miRNAs targets, Peaks from ChIP-seq, ATAC, et al.)  The first three columns should be chromosome, start, end. Example: 1   567821 573421  EEF1D

        -g    GWAS summary statistics. The first two columns should be chromosome and position, the last column shoud be effect values,such as t value or beta value. Example: 1  123089 rs0011345  0.00045  -1.4625

        -e    extended range (KB) for genome feature regions. For example, the -e 100 means  genome feature regions should also include their unstream/downstream 100kb region. Default -e is 0.

        -n    repeat n times for the permutation test. Default -n is 1000.

        -o    output file. Defualt the output file name is \"SumGSE_permutation.out\".

        -s    specific regions (Optional). If assuming -b, the permutation SNPs will be limited in these specific regions. 

";


foreach my $i (0 ..scalar(@ARGV)-1) {
  if($ARGV[$i] eq '-i') {
    $inputfile = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-g') {
    $GWAS = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-s') {
    $specific = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-n'){
	$times = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-e'){
        $extention = $ARGV[++$i]*1000;
  }elsif($ARGV[$i] eq '-o'){
        $outfile = $ARGV[++$i];
   }elsif($ARGV[$i]=~ m/\-/){
    die "Cannot recognize the parameter! $ARGV[$i].\n$usage";
   }
}

$extention ||= 0;
$times ||= 1000;
$outfile ||= "SumGSE_permutation.out";
test_first_argument();

my %hash_chr;
my %gene_pos;
my $link='_';
    open (T1,"$inputfile") or die "Cannot open the input file! $inputfile.\n$usage";
    while($line=<T1>){
    chomp $line;
    @terms =split (/\t/ ,$line);
    my($chr,$start,$end) = @terms[0,1,2];
    $gene_pos{$chr.$link.$start}=$end;  
    if(exists $hash_chr{$chr}){
        my $ref = $hash_chr{$chr};
        push(@{$$ref[1]}, $start);
        push(@{$$ref[2]}, $end);
    }else{
        my @starts = ($start);
        my @ends = ($end);
        my $data = [$chr,\@starts,\@ends];
        $hash_chr{$chr} = $data;
    }
    }
    close(T1);

my %SNP_chr;
my %SNP_pos;
my %SNP_t;
my %rank;
my $within=0;
my $SNP_num=0;
my $value=0;
my @within_snps;
my $count=0;
my $sum_per=0;   
if($specific ne ''){
my %gene_pos2;
my %hash_chr2;
        open (T3,"$specific") or die "Cannot open specific regions file(-s) $specific.\n$usage";
        while($line2=<T3>){
            chomp $line2;
            @terms2 =split (/\t/ ,$line2);
            my($chr2,$start2,$end2) = @terms2[0,1,2];
            $gene_pos2{$chr2.$link.$start2}=$end2;
            if(exists $hash_chr2{$chr2}){
            my $ref2 = $hash_chr2{$chr2};
            push(@{$$ref2[1]}, $start2);
            push(@{$$ref2[2]}, $end2);
            }else{
            my @starts2 = ($start2);
            my @ends2 = ($end2);
            my $data2 = [$chr2,\@starts2,\@ends2];
            $hash_chr2{$chr2} = $data2;
            }
        }
        close(T3);

        open (T2,"$GWAS") or die "Cannot open $GWAS.\n$usage";
        while($line=<T2>){
            chomp $line;
            $SNP_num++;
            @A =split (/\t/ ,$line);
            my($pos_chr,$pos,$t) = @A[1,2,-1];
            $SNP_t{$pos_chr.$link.$pos}=$t;
            if(exists $hash_SNP{$pos_chr}){
            my $snp = $hash_SNP{$pos_chr};
            push(@{$$snp[1]}, $pos);
            push(@{$$snp[2]}, $t);
            }else{
            my @poss = ($pos);
            my @ts = ($t);
            my $data3 = [$pos_chr,\@poss,\@ts];
            $hash_SNP{$pos_chr} = $data3;
            }
        }
        close(T2);

my %SNP_t1;
my $value=0;
my $within_allgene=0;
        foreach my $key2 (keys %hash_chr2) {
            my $ref2 = $hash_chr2{$key2};
            my @starts2 = @{$$ref2[1]};
            my @ends2 = @{$$ref2[2]};
            if (exists $hash_SNP{$key2}){
                my $ref3 = $hash_SNP{$key2};
                my @poss3 = @{$$ref3[1]};
                my @tvs3 = @{$$ref3[2]};
                foreach my $pos3 (@poss3){
                    foreach my $start2 (@starts2){
                        if($pos3>=$start2-$extention && $pos3<=$gene_pos2{$key2.$link.$start2}+$extention){
                            if(exists $rank{$key2.$link.$pos3}){
                                next;}else{
                                $within_allgene++;
                                $rank{$key2.$link.$pos3}=$within_allgene;
                                }
                            }
                    }
                }
             }
        }
open(my $ot, '>', $outfile) or die "Could not open file '$outfile' $!";
print "The number of SNPs within specific regions"."\t";
print "$within_allgene\n";

my $value1=0;
my %hash2;
 foreach my $key (keys %hash_chr) {
        my $ref = $hash_chr{$key};
        my @starts = @{$$ref[1]};
        my @ends = @{$$ref[2]};
                if (exists $hash_SNP{$key}){
                my $newref = $hash_SNP{$key};
                my @poss = @{$$newref[1]};
                my @tvs = @{$$newref[2]};
                        foreach my $posnew (@poss){
                                foreach my $startnew (@starts){
                                        if($posnew>=$startnew-$extention && $posnew<=$gene_pos{$key.$link.$startnew}+$extention){
                                         if(exists $hash2{$key.$link.$posnew}){
                    next;}else{ 
                    $within++;
                    $hash2{$key.$link.$posnew}="";
                                        $value1=$value1+$SNP_t{$key.$link.$posnew}**2;
                                        push (@within_snps, $rank{$key.$link.$posnew});
                                        }
                                }
                        }
                }
        }
}
print  "All SNP number in GWAS summary"."\t";
print  "$SNP_num\n";
print  "SNPs within feature regions"."\t";
print  "$within\n";
print $ot "ID"."\t";
print $ot "sum(tvalue*2)"."\n";
print $ot "Feature"."\t";
print  $ot "$value1\n";  
  %rank_reverse=reverse (%rank);
my $count=0;
my $sum_per=0;
for( my $a = 1; $a <=$times; $a = $a + 1 ){
    my $num=int(rand($within_allgene))+1;
    my @new_rank;
    foreach my $perm (@within_snps){
        if($perm<=$within_allgene-$num){
            $new_rank=$perm+$num;
        }else{
            $new_rank=$perm+$num-$within_allgene;
        }
           push(@new_rank,$new_rank);
    }
#print join(",",@new_rank);
my $zzz;
my $per_value=0;
foreach my $new_id (@new_rank){
    $zzz=$rank_reverse{$new_id};
#   print "$zzz\n";
    $per_value=$per_value+$SNP_t{$zzz}**2;
    }
print $ot "Permutation_";  
print $ot "$a\t";
print $ot "$per_value\n";
if($per_value>$value1){
    $count++;
    }
 $sum_per=$sum_per+$per_value;
}
my $pvalue=($count+1)/($times+1);
my $fold=($value1*$times)/$sum_per;
print "P-value"."\t";
print "$pvalue\n";
print "Fold"."\t";
print "$fold\n";

print $ot "P-value"."\t";
print $ot "$pvalue\n";
print $ot "Fold"."\t";
print $ot "$fold\n";
}






else{
my %SNP_chr;
my %SNP_pos;
my %SNP_t;
my %rank;
my $within=0;
my $SNP_num=0;
my $value=0;
my @within_snps;
        open (T2,"$GWAS") or die "Cannot open $GWAS.\n$usage";
        while($line=<T2>){
        chomp $line;
    $SNP_num++;
        @A =split (/\t/ ,$line);
    my($pos_chr,$pos,$t) = @A[1,2,-1];
    $SNP_t{$pos_chr.$link.$pos}=$t;
     $rank{$pos_chr.$link.$pos}=$SNP_num;
    if(exists $hash_SNP{$pos_chr}){
        my $snp = $hash_SNP{$pos_chr};
        push(@{$$snp[1]}, $pos);
        push(@{$$snp[2]}, $t);
        }else{
        my @poss = ($pos);
        my @ts = ($t);
        my $data3 = [$pos_chr,\@poss,\@ts];
        $hash_SNP{$pos_chr} = $data3;
        }
    }
close(T2);

#my $extention="$ARGV[2]";
        #$SNP_t{$A[1]}=$A[2];
my $value1=0;
my %hash2;
 foreach my $key (keys %hash_chr) {
        my $ref = $hash_chr{$key};
        my @starts = @{$$ref[1]};
        my @ends = @{$$ref[2]};
                if (exists $hash_SNP{$key}){
                my $newref = $hash_SNP{$key};
                my @poss = @{$$newref[1]};
                my @tvs = @{$$newref[2]};
                        foreach my $posnew (@poss){
                                foreach my $startnew (@starts){
                                        if($posnew>=$startnew-$extention && $posnew<=$gene_pos{$key.$link.$startnew}+$extention){
                                         if(exists $hash2{$key.$link.$posnew}){
                    next;}else{ 
                    $within++;
                    $hash2{$key.$link.$posnew}="";
                                        $value1=$value1+$SNP_t{$key.$link.$posnew}**2;
                                        push (@within_snps, $rank{$key.$link.$posnew});
                                        }
                                }
                        }
                }
        }
}
open(my $ot, '>', $outfile) or die "Could not open file '$outfile' $!";
print "All SNP number in GWAS summary"."\t";
print "$SNP_num\n";
print "SNPs within feature regions"."\t";
print "$within\n";
print $ot "ID"."\t";
print $ot "sum(tvalue*2)"."\n";
print $ot "Feature"."\t";
print $ot "$value1\n";

%rank_reverse=reverse (%rank);
my $count=0;
my $sum_per=0;
for( my $a = 1; $a <=$times; $a = $a + 1 ){
    my $num=int(rand($SNP_num))+1;
    my @new_rank;
    foreach my $perm (@within_snps){
        if($perm<=$within-$num){
            $new_rank=$perm+$num;
        }else{
            $new_rank=$perm+$num-$within;
        }
           push(@new_rank,$new_rank);
    }
my $zzz;
my $per_value=0;
foreach my $new_id (@new_rank){
    $zzz=$rank_reverse{$new_id};
    $per_value=$per_value+$SNP_t{$zzz}**2;
    }
print $ot "Permutation_";  
print $ot "$a\t";
print $ot "$per_value\n";
if($per_value>$value1){
    $count++;
    }
 $sum_per=$sum_per+$per_value;
}
my $pvalue=($count+1)/($times+1);
my $fold=($value1*$times)/$sum_per;
print $ot "P-value"."\t";
print $ot "$pvalue\n";
print $ot "Fold"."\t";
print $ot "$fold\n";

print "P-value"."\t";
print "$pvalue\n";
print "Fold"."\t";
print "$fold\n";
}

sub test_first_argument{
    if(not $ARGV[0]){die $usage;}
    if($ARGV[0] eq '-h' or $ARGV[0] eq '--help'){die $usage;}
    return;
   # if(not $extention){$extention=0;}
}

close($ot);
