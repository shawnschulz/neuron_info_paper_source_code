#!/bin/bash
# 

super_cluster=$1

mkdir -p ${super_cluster}_fastas
mv $super_cluster*_cell_fasta.fa ${super_cluster}_fastas/
