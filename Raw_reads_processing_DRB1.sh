#####Asses sequence quality (FASTQC) 
fastqc 02001_S1_L001_R1_001.fastq.gz 
fastqc 02001_S1_L001_R2_001.fastq.gz

####Concatenate reads
cat 02001_S1_L001_R1_001.fastq.gz 02001_S1_L001_R2_001.fastq.gz > 02001_S1_L001_001.fastq.gz
gzip -d 02001_S1_L001_001.fastq.gz

##Filter reads by quality with FASTX (FASTX-TOOLKIT)
fastq_quality_filter -q 20 -p 99 -i 02001_S1_L001_001.fastq -o 02001_filtered.fastq

##Demultiplex DRB1 and cut 5' primer (CUTADAPT)
cutadapt --overlap 15 -g DRB1_f=GTCCCCACAG -g DRB1_r=GCGCGGAGTCTCGGCAGGG 02001_filtered.fastq -o 02001_filtered_5_prime_trimmed-{name}.fastq

##Cut 3' primer (CUTADAPT)
cutadapt --overlap 15 -a DRB1_f=GTGAGCGCGGGGCTGGGCGGG 02001_filtered_5_prime_trimmed-DRB1_f.fastq -o 02001_filtered_5_prime_trimmed-DRB1_f{name}.fastq

#####Cut introns (CUTADAPT)
cutadapt --overlap 15 -e 0.21 -a CCGCAGAG...GGTGAGCGCAGGCCGCCCTCCGCGGGGCCCACCCTC 02001_merged_filtered2_both_trimmed_remerged-DRB1.fastq -o 02001_merged_filtered2_both_trimmed_remerged-DRB1_no_intron.fastq

####Final filter (FASTX-TOOLKIT)
fastq_quality_filter -q 30 -p 95 -i 02001_merged_filtered2_both_trimmed_remerged-DRB1_no_intron.fastq -o 02001_DRB1_final.fastq
