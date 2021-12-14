#!/usr/bin/env bash

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.agp.gz | gunzip - > hg38.agp
awk 'BEGIN {OFS="\t";} {component_type=$5; chr=$1; start=$2-1; end=$3;} {if((component_type=="N" || component_type=="U") && chr!~/_alt/ && chr!~/chrUn_/ && chr!~/_random/) print chr, start, end}' hg38.agp | sort -k1,1 -k2,2n > hg38.gaps.bed
awk 'BEGIN {OFS="\t";} {chr=$1; size=$2;} {if(chr!~/_alt/ && chr!~/chrUn_/ && chr!~/_random/) print chr, size}' hg38.chrom.sizes |\
sort -k1,1 >\
hg38-main.chrom.sizes

# make windows
bedtools makewindows -g hg38-main.chrom.sizes -w 20000 > hg38.20kb.windows.bed
# remove windows with any gap and any target
awk 'BEGIN {OFS="\t";} {chr=$1; start=$2-10000; end=$3+10000;} { print chr, start, end}' Agilent_V5_targets.bed > Agilent_V5_targets.up10kb.bed
bedtools intersect -v -a hg38.20kb.windows.bed -b Agilent_V5_targets.up10kb.bed hg38.gaps.bed > hg38.20kb.windows-offtarget.bed
