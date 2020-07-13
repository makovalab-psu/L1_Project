####################################################################################################
#Analyze genome-wide distribution of polyA/T sequences##
####################################################################################################

#####Examined the number of the de novo L1 element reads and their 1-kb flanks that overlap with the genomic polyA/T sequences #####

cd ~/Google_Drive/L1_Project/Analysis/Additional_analysis/PolyAT_analysis

#Download hg19 genome sequence (hg.fa.gz)
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
#Upload sequences to Galaxy and extract polyA locations using "Nucleotide subsequence search" function: 
#[Nucleotide subsequence search providing regions in BED format (Galaxy Version 0.2)]
#Shared history
#https://usegalaxy.org/u/cdbruce/h/ployahg19genome
#Download .bed format output containing genomic coordinates of all polyA seqeuences (>15 bp)
#Remove non-standard chromosomes and Y
cat Galaxy5-[polyA_15_hg19].bed | awk  '$1 !~ /_/ && $1 !~ /M/ && $1 !~ /Y/' | cut -f1,2,3,5 | sortBed -i stdin > polyA_hg19_clean_sorted.bed

#Extract insertion site for each de novo L1 dataset (by strand)
cat L1denovo_BWA_17037_reads.bed | cut -f1,2,3,4 | awk '{ if ($4 == "+") print $1"\t"$2"\t"$2+1"\t"$4;  if ($4 == "-") print $1"\t"$3"\t"$3+1"\t"$4}' > L1denovo_BWA_17037_insertions.bed

#Overlap bettwem de novo L1s (all 17037) and polyA
cat L1denovo_BWA_17037_reads.bed | cut -f1,2,3,4 | bedtools coverage -a stdin -b polyA_hg19_clean_sorted.bed > cov_denovo_polyA.bed
wc -l L1denovo_BWA_17037_reads.bed
#17037
cat cov_denovo_polyA.bed | awk '{ if ($8 != 0) print $0}' | wc -l
#7

#Overlap bettwem de novo L1 1-kb flanks and polyA
cat L1denovo_BWA_17037_reads.bed | cut -f1,2,3,4 | awk '{ print $1"\t"($2-500)"\t"($3+500)"\t"$4}' | sortBed > L1denovo_BWA_17037_1kb_flanks.bed
wc -l L1denovo_BWA_17037_1kb_flanks.bed
#17037
bedtools coverage -a L1denovo_BWA_17037_1kb_flanks.bed -b polyA_hg19_clean_sorted.bed | awk '{ if ($8 != 0) print $0}' | wc -l
#207
