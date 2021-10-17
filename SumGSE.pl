#!/usr/bin/perl
###########################################################################
## SumGSE.pl --- Process genome region bed,GWAS summary
############################################################################
## Author: Wentao Cai <wtaocai@gmail.com>
## updated: 17 Oct 2021
## Version: 0.02
my $usage = "\nUsage:\n\t perl SumGSE.pl -a [genome_region.bed] -b [GWAS_summaries.txt] -g [specific.regions] -e [genome_region_extention] -n [permutation_times]\n\n";
foreach my $i (0 ..scalar(@ARGV)-1) {
  if($ARGV[$i] eq '-a') {
    $T1 = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-b') {
    $T2 = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-g') {
    $T3 = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-n'){
	$times = $ARGV[++$i];
 }elsif($ARGV[$i] eq '-e'){
        $extention = $ARGV[++$i];
  }

if(@ARGV ==0) {
    die $usage;
}
}


my %hash_chr;
my %gene_pos;
my $link='_';
        open (T1,"$T1") or die "Cannot open $T1.\n$usage";
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

my %gene_pos2;
my %hash_chr2;
        open (T3,"$T3") or die "Cannot open $T3.\n$usage";
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

my %SNP_chr;
my %SNP_pos;
my %SNP_t;
my %rank;
my $within=0;
my $SNP_num=0;
my $value=0;
my @within_snps;
        open (T2,"$T2") or die "Cannot open $T2.\n$usage";
        while($line=<T2>){
        chomp $line;
	$SNP_num++;
        @A =split (/\t/ ,$line);
	my($pos_chr,$pos,$t) = @A[1,2,-1];
       if($A[1]==14){
       if($A[2] > 795351 && $A[2]<2804562){
       next;
       }
       }
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

#my $extention="$ARGV[2]";
        #$SNP_t{$A[1]}=$A[2];
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
                                        #$value=$value+$SNP_t{$key2.$link.$pos3}**2;
					$rank{$key2.$link.$pos3}=$within_allgene;
         				#$SNP_t1{$key2.$link.$pos3}=$SNP_t{$key2.$link.$pos3}                             
         #                               push (@within_snps, $rank{$key.$link.$pos1});
                                }
				}
                        }
                }
        }
}
print "The number of SNP within genes"."\t";
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
print "SNP_number"."\t";
print "$SNP_num\n";
print "SNP_number_within_feature_region"."\t";
print "$within\n";
print "sum(tvalue*2)"."\t";
print "$value1\n";	

#my $times= "$ARGV[3]";
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
#	print "$zzz\n";
	$per_value=$per_value+$SNP_t{$zzz}**2;
	}
print "$a\t";
print "$per_value\n";
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


