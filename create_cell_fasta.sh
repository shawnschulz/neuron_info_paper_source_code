#!/bin/bash
#
ensembl_id=$1
file_name=$2
number_reads=$3

cd ~/Programming/neuron_rna_info_paper/data

for i in $(seq 0 $number_reads); do
	pueue add "echo ${ensembl_id}_${i} >> $file_name; grep -A1 $ensembl_id ../data/ensembl_annotated_transcriptome.fa | grep -v $ensembl_id >> $file_name"
done
