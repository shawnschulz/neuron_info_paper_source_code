#!/bin/bash
#
for i in sc*; do
    name=${i%.h5ad}
    command="python3 /home/shawn/Programming/neuron_rna_info_paper/src/get_counts.py --input $i --name $name"
    pueue add "$command"
done
