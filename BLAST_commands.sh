#!/bin/bash

#Load BLAST
ml bio/BLAST+/2.11.0-gompi-2020b
ml lang/Anaconda3/2021.05
source activate myenv

#BLAST MHC exons to each Balaenoptera_Acutorostrata genome and outputs a file with 12 columns (query ID, scaffold of subject, percentage identical matches, alignment length, number of mismatches, number of gap openings, query start, query end, subject start, subject end, evalue, bitscore)
blastn -query /work/celeminamaro/celeminamaro/celeminamaro/BLAST_MHC_Cetaceans_Niko/query.txt -db Balaenoptera_Acutorostrata.db -out Balaenoptera_Acutorostrata.bed -outfmt "6 qseqid sseqid sstart send sframe"

#Edit start/end position of the sequences in 3' to 5' direction to be in the correct order
awk '{if ($5=="-1") {tmp=$3; $3=$4; $4=tmp;} print $0}' Balaenoptera_Acutorostrata.bed > Balaenoptera_Acutorostrata_corrected.bed

#Prepare .bed file for bedtools
cat Balaenoptera_Acutorostrata_corrected.bed | awk '{print $2 "\t" $3 "\t" $4}' > Balaenoptera_Acutorostrata_corrected_coords.bed
awk '$3+=100' Balaenoptera_Acutorostrata_corrected_coords.bed > Balaenoptera_Acutorostrata_corrected_coords2.bed
awk '$2-=100' Balaenoptera_Acutorostrata_corrected_coords2.bed > Balaenoptera_Acutorostrata_corrected_coords3.bed
awk -v OFS="\t" '$1=$1' Balaenoptera_Acutorostrata_corrected_coords3.bed > Balaenoptera_Acutorostrata_corrected_coords4.bed

#Extract blast hits from the genome as fasta sequences
bedtools getfasta -fi GCA_949987535.1_mBalAcu1.1_genomic.fna -bed Balaenoptera_Acutorostrata_corrected_coords4.bed -fo Balaenoptera_Acutorostrata_MHCII_exon2.fasta





