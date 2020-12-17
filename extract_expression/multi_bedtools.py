source activate python36
bedtools multicov -bams \
0619_B--TACAGC_1_tripped_sort.bam \
0619_A--ACTTGA_1_tripped_sort.bam \
3008_C--ACTGAT_1_tripped_sort.bam \
0619_D--GGCTAC_1_tripped_sort.bam \
3008_A--GTCCGC_1_tripped_sort.bam \
3678_A--GTGAAA_1_tripped_sort.bam \
3678_D--AGTCAA_1_tripped_sort.bam \
3678_D--CACTCA_1_tripped_sort.bam \
3678_B--CGTACG_1_tripped_sort.bam \
3008_A--ATTCCT_1_tripped_sort.bam \
0619_B--GTGGCC_1_tripped_sort.bam \
3008_C--TCCCGA_1_tripped_sort.bam \
3678_C--TAGCTT_1_tripped_sort.bam \
3678_C--CAACTA_1_tripped_sort.bam \
3008_D--CACGAT_1_tripped_sort.bam \
3678_B--TCATTC_1_tripped_sort.bam \
3008_B--TATAAT_1_tripped_sort.bam \
0619_D--CACCGG_1_tripped_sort.bam \
3008_B--GTTTCG_1_tripped_sort.bam \
3678_A--TAATCG_1_tripped_sort.bam \
0619_A--CCGTCC_1_tripped_sort.bam \
0619_C--GAGTGG_1_tripped_sort.bam \
3008_D--CTTGTA_1_tripped_sort.bam \
0619_C--CAGATC_1_tripped_sort.bam \
-bed ./hsa.gff3 > result.txt