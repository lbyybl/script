#!/bin/bash

# Boyuan_Li

# 2019/11/14

bedtools intersect -a dox05hfirst_high_deseq2.txt -b R05h1first_high_deseq2.txt -wa -wb | awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' > overlap/jj.bed
bedtools intersect -a dox05hfirst_high_deseq2.txt -b R05h1unchange_deseq2.txt -wa -wb | awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' > overlap/jp.bed
bedtools intersect -a dox05hfirst_high_deseq2.txt -b R05h1second_high_deseq2.txt -wa -wb | awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' > overlap/jz.bed
bedtools intersect -a dox05hunchange_deseq2.txt -b R05h1first_high_deseq2.txt -wa -wb | awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' > overlap/pj.bed
bedtools intersect -a dox05hunchange_deseq2.txt -b R05h1unchange_deseq2.txt -wa -wb | awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' > overlap/pp.bed
bedtools intersect -a dox05hunchange_deseq2.txt -b R05h1second_high_deseq2.txt -wa -wb | awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' > overlap/pz.bed
bedtools intersect -a dox05hsecond_high_deseq2.txt -b R05h1first_high_deseq2.txt -wa -wb | awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' > overlap/zj.bed
bedtools intersect -a dox05hsecond_high_deseq2.txt -b R05h1unchange_deseq2.txt -wa -wb | awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' > overlap/zp.bed
bedtools intersect -a dox05hsecond_high_deseq2.txt -b R05h1second_high_deseq2.txt -wa -wb | awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' > overlap/zz.bed
ls overlap/*.bed | parallel "bedtools merge -i <(sort -k 1,1 -k 2n,2 -k 3n,3 {}) | sort -k 1,1 -k 2n,2 -k 3n,3 -o {}"
/home/boyuanli/bashscript/bin/ruvb2/chip-seq/rm_ovalap.sh -d overlap

