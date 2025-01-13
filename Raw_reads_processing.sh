#####Asses sequence quality (FASTQC) 
fastqc 02001_S1_L001_R1_001.fastq.gz 
fastqc 02001_S1_L001_R2_001.fastq.gz

#####Merge reads (FLASH2)
flash2 02001_S1_L001_R1_001.fastq.gz 02001_S1_L001_R2_001.fastq.gz -o 02001

#####First quality filter (FASTX-TOOLKIT)
fastq_quality_filter -q 20 -p 90 -i 02001.extendedFrags.fastq -o 02001_merged_filtered.fastq

#####Demultiplex and cut 5' primer (CUTADAPT)
cutadapt --overlap 15 -g DRA_f=CTCCTGGTTCCCACCCTAAT -g DRA_r=TGGCATGAGAATTTGGAGCT 02001_merged_filtered2.fastq -o 02001_merged_filtered2_5_prime_trimmed-{name}.fastq
cutadapt --overlap 15 -g DRB2_f=CTGCACCGTGAAGCACTCCA -g DRB2_r=GTGTCTCCCGCAGAGAATTT 02001_merged_filtered2.fastq -o 02001_merged_filtered2_5_prime_trimmed-{name}.fastq
cutadapt --overlap 15 -g DQB_f=CTGGCTGAGCGGCGGTGTCTC -g DQB_r=GCGCGGAGTCTCGGCAGGG 02001_merged_filtered2.fastq -o 02001_merged_filtered2_5_prime_trimmed-{name}.fastq
cutadapt --overlap 15 -g DQA_f=CTGCCCCTCACCTTCACTTA -g DQA_r=GGGAACAAGAGAGTGAGGCC 02001_merged_filtered2.fastq -o 02001_merged_filtered2_5_prime_trimmed-{name}.fastq
cutadapt --overlap 15 -g DOA_f=GGCACAGCTCAGAGACTCAA -g DOA_r=CAGGAGTCAACAGAGGCCTG 02001_merged_filtered2.fastq -o 02001_merged_filtered2_5_prime_trimmed-{name}.fastq
cutadapt --overlap 15 -g DMA_f=GCAGGCCTTGCTGTTTTCTA -g DMA_r=TTAAGAGGGCAGCCCACCCC 02001_merged_filtered2.fastq -o 02001_merged_filtered2_5_prime_trimmed-{name}.fastq
cutadapt --overlap 15 -g DMB_f=GACCTATCCCCCTTCTCTGG -g DMB_r=CCCTCTGCTCCCTTCTCGCC 02001_merged_filtered2.fastq -o 02001_merged_filtered2_5_prime_trimmed-{name}.fastq
cutadapt --overlap 15 -g DOB_f=GTTTCCTAATTGGGGCTGGT -g DOB_r=GTCCCCACCTTACACCTCAC 02001_merged_filtered2.fastq -o 02001_merged_filtered2_5_prime_trimmed-{name}.fastq

#####Cut 3' primer on reverse reads (CUTADAPT)
cutadapt --overlap 15 -a DRA_r=AGCTCCAAATTCTCATGCCA 02001_merged_filtered2_5_prime_trimmed-DRA_f.fastq -o 02001_merged_filtered2_both_trimmed-DRA_f-{name}.fastq
cutadapt --overlap 15 -a DRB2_r=AAATTCTCTGCGGGAGACAC 02001_merged_filtered2_5_prime_trimmed-DRB2_f.fastq -o 02001_merged_filtered2_both_trimmed-DRB2_f-{name}.fastq
cutadapt --overlap 15 -a DQB_r=CCCTGCCGAGACTCCGCGC 02001_merged_filtered2_5_prime_trimmed-DQB_f.fastq -o 02001_merged_filtered2_both_trimmed-DQB_f-{name}.fastq
cutadapt --overlap 15 -a DQA_r=GGCCTCACTCTCTTGTTCCC 02001_merged_filtered2_5_prime_trimmed-DQA_f.fastq -o 02001_merged_filtered2_both_trimmed-DQA_f-{name}.fastq
cutadapt --overlap 15 -a DOA_r=CAGGCCTCTGTTGACTCCTG 02001_merged_filtered2_5_prime_trimmed-DOA_f.fastq -o 02001_merged_filtered2_both_trimmed-DOA_f-{name}.fastq
cutadapt --overlap 15 -a DMA_r=GGGGTGGGCTGCCCTCTTAA 02001_merged_filtered2_5_prime_trimmed-DMA_f.fastq -o 02001_merged_filtered2_both_trimmed-DMA_f-{name}.fastq
cutadapt --overlap 15 -a DMB_r=GGCGAGAAGGGAGCAGAGGG 02001_merged_filtered2_5_prime_trimmed-DMB_f.fastq -o 02001_merged_filtered2_both_trimmed-DMB_f-{name}.fastq
cutadapt --overlap 15 -a DOB_r=GTGAGGTGTAAGGTGGGGAC 02001_merged_filtered2_5_prime_trimmed-DOB_f.fastq -o 02001_merged_filtered2_both_trimmed-DOB_f-{name}.fastq

#####Cut 3' primer on forward reads (CUTADAPT)
cutadapt --overlap 15 -a DRA_f=ATTAGGGTGGGAACCAGGAG 02001_merged_filtered2_5_prime_trimmed-DRA_r.fastq -o 02001_merged_filtered2_both_trimmed-DRA_r-{name}.fastq
cutadapt --overlap 15 -a DRB2_f=TGGAGTGCTTCACGGTGCAG 02001_merged_filtered2_5_prime_trimmed-DRB2_r.fastq -o 02001_merged_filtered2_both_trimmed-DRB2_r-{name}.fastq
cutadapt --overlap 15 -a DQB_f=GAGACACCGCCGCTCAGCCAG 02001_merged_filtered2_5_prime_trimmed-DQB_r.fastq -o 02001_merged_filtered2_both_trimmed-DQB_r-{name}.fastq
cutadapt --overlap 15 -a DQA_f=TAAGTGAAGGTGAGGGGCAG 02001_merged_filtered2_5_prime_trimmed-DQA_r.fastq -o 02001_merged_filtered2_both_trimmed-DQA_r-{name}.fastq
cutadapt --overlap 15 -a DOA_f=TTGAGTCTCTGAGCTGTGCC 02001_merged_filtered2_5_prime_trimmed-DOA_r.fastq -o 02001_merged_filtered2_both_trimmed-DOA_r-{name}.fastq
cutadapt --overlap 15 -a DMA_f=TAGAAAACAGCAAGGCCTGC 02001_merged_filtered2_5_prime_trimmed-DMA_r.fastq -o 02001_merged_filtered2_both_trimmed-DMA_r-{name}.fastq
cutadapt --overlap 15 -a DMB_f=CCAGAGAAGGGGGATAGGTC 02001_merged_filtered2_5_prime_trimmed-DMB_r.fastq -o 02001_merged_filtered2_both_trimmed-DMB_r-{name}.fastq
cutadapt --overlap 15 -a DOB_f=ACCAGCCCCAATTAGGAAAC 02001_merged_filtered2_5_prime_trimmed-DOB_r.fastq -o 02001_merged_filtered2_both_trimmed-DOB_r-{name}.fastq

#####Reverse complement reverse reads. Previously checked that were not in ORF and that forward reads were in ORF (SEQTK) 
seqtk seq -r 02001_merged_filtered2_both_trimmed-DRA_r-DRA_f.fastq > 02001_merged_filtered2_both_trimmed-DRA_r_RC.fastq
seqtk seq -r 02001_merged_filtered2_both_trimmed-DRB2_r-DRB2_f.fastq > 02001_merged_filtered2_both_trimmed-DRB2_r_RC.fastq
seqtk seq -r 02001_merged_filtered2_both_trimmed-DQB_r-DQB_f.fastq > 02001_merged_filtered2_both_trimmed-DQB_r_RC.fastq
seqtk seq -r 02001_merged_filtered2_both_trimmed-DQA_r-DQA_f.fastq > 02001_merged_filtered2_both_trimmed-DQA_r_RC.fastq
seqtk seq -r 02001_merged_filtered2_both_trimmed-DOA_r-DOA_f.fastq > 02001_merged_filtered2_both_trimmed-DOA_r_RC.fastq
seqtk seq -r 02001_merged_filtered2_both_trimmed-DMA_r-DMA_f.fastq > 02001_merged_filtered2_both_trimmed-DMA_r_RC.fastq
seqtk seq -r 02001_merged_filtered2_both_trimmed-DMB_r-DMB_f.fastq > 02001_merged_filtered2_both_trimmed-DMB_r_RC.fastq
seqtk seq -r 02001_merged_filtered2_both_trimmed-DOB_r-DOB_f.fastq > 02001_merged_filtered2_both_trimmed-DOB_r_RC.fastq

#####Merge forward and reverse (in reverse complement) reads
cat 02001_merged_filtered2_both_trimmed-DRA_f-DRA_r.fastq 02001_merged_filtered2_both_trimmed-DRA_r_RC.fastq > 02001_merged_filtered2_both_trimmed_remerged-DRA.fastq
cat 02001_merged_filtered2_both_trimmed-DRB2_f-DQB_r.fastq 02001_merged_filtered2_both_trimmed-DRB2_r_RC.fastq > 02001_merged_filtered2_both_trimmed_remerged-DRB2.fastq
cat 02001_merged_filtered2_both_trimmed-DQB_f-DQB_r.fastq 02001_merged_filtered2_both_trimmed-DQB_r_RC.fastq > 02001_merged_filtered2_both_trimmed_remerged-DQB.fastq
cat 02001_merged_filtered2_both_trimmed-DQA_f-DQA_r.fastq 02001_merged_filtered2_both_trimmed-DQA_r_RC.fastq > 02001_merged_filtered2_both_trimmed_remerged-DQA.fastq
cat 02001_merged_filtered2_both_trimmed-DOA_f-DOA_r.fastq 02001_merged_filtered2_both_trimmed-DOA_r_RC.fastq > 02001_merged_filtered2_both_trimmed_remerged-DOA.fastq
cat 02001_merged_filtered2_both_trimmed-DMA_f-DMA_r.fastq 02001_merged_filtered2_both_trimmed-DMA_r_RC.fastq > 02001_merged_filtered2_both_trimmed_remerged-DMA.fastq
cat 02001_merged_filtered2_both_trimmed-DMB_f-DMB_r.fastq 02001_merged_filtered2_both_trimmed-DMB_r_RC.fastq > 02001_merged_filtered2_both_trimmed_remerged-DMB.fastq
cat 02001_merged_filtered2_both_trimmed-DOB_f-DOB_r.fastq 02001_merged_filtered2_both_trimmed-DOB_r_RC.fastq > 02001_merged_filtered2_both_trimmed_remerged-DOB.fastq

####Final filter (FASTX-TOOLKIT)
fastq_quality_filter -q 30 -p 95 -i 02001_merged_filtered2_both_trimmed_remerged-DRA.fastq -o 02001_DRA_final.fastq
fastq_quality_filter -q 30 -p 95 -i 02001_merged_filtered2_both_trimmed_remerged-DRB2.fastq -o 02001_DRB2_final.fastq
fastq_quality_filter -q 30 -p 95 -i 02001_merged_filtered2_both_trimmed_remerged-DQB.fastq -o 02001_DQB_final.fastq
fastq_quality_filter -q 30 -p 95 -i 02001_merged_filtered2_both_trimmed_remerged-DQA.fastq -o 02001_DQA_final.fastq
fastq_quality_filter -q 30 -p 95 -i 02001_merged_filtered2_both_trimmed_remerged-DOA.fastq -o 02001_DOA_final.fastq
fastq_quality_filter -q 30 -p 95 -i 02001_merged_filtered2_both_trimmed_remerged-DMA.fastq -o 02001_DMA_final.fastq
fastq_quality_filter -q 30 -p 95 -i 02001_merged_filtered2_both_trimmed_remerged-DMB.fastq -o 02001_DMB_final.fastq
fastq_quality_filter -q 30 -p 95 -i 02001_merged_filtered2_both_trimmed_remerged-DOB.fastq -o 02001_DOB_final.fastq
