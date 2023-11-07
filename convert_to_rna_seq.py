#!/usr/bin/env python
import scanpy as sc 
import os
import subprocess

os.chdir("../data")

os.getcwd()

def write_to_file(string, fn):
    with open(fn, 'w') as file:
        file.write(string)
    file.close()
    return(1)

def multiply_lines_with_substring(multi_line_string, substring, n):
    #this is waaaaay too slow
    lines = multi_line_string.split('\n')
    result_lines = []
    for line in lines:
        if substring in line:
            result_lines.extend([line] * n)
        else:
            result_lines.append(line)
    result_string = '\n'.join(result_lines)
    return result_string

#test to look how much memory is being used by write string
def utf8len(s):
    return len(s.encode('utf-8'))

adata = sc.read_h5ad("sc_splatter.h5ad")

#testing on just 10 cells for now
adata = adata[adata.obs_names[0:1]]

ref_fasta_path = "ensembl_annotated_transcriptome.fa"

super_cluster_term = "splatter_one_cell"
dir_name = super_cluster_term + "_fastas"
subprocess.call(["mkdir", "-p", dir_name])

out_fasta_fn =  super_cluster_term + '.fa' 

string_to_write = "" #this string is gonna be REALLY long...

cell_iteration_count=0
gene_index_chosen_count=0

subprocess.run(["echo", "", ">", out_fasta_fn])

#should spin up threads and just run all of the cell_indexes in a parallelized queue
#lets refactor to make this happen
for cell_index in range(0, len(adata.obs_names)): 
    print("DEBUG: cell iterations is: {}".format(cell_iteration_count))
    cell_iteration_count+=1
    print("The current cell index is {}".format(cell_index))
    cell_name = adata.obs_names[cell_index]

    rg_args = ""
    for gene_index in range(0, len(adata.var_names)):
        print("The current gene_index is {}".format(gene_index))
        read_count = int(adata.X[cell_index, gene_index])
        ensembl_id = "gene:" + adata.var_names[gene_index]
        #should cross check the p-value of expression of the gene, hopefully it's already in the adata
        #we can use that to figure out how likely the gene is to actually be differentially expressed 
        if read_count > 0:
            gene_index_chosen_count += 1
            print("DEBUG: gene index chosen count is: {}".format(gene_index_chosen_count))
            rg_args = rg_args + ensembl_id + '\n'  
        else:
            print("no reads associated with this transcript")
            
    print("DEBUG: printing rg_args: {}".format(rg_args))
    rg_args_fn = "temp_rg_args.txt"
    write_to_file(rg_args, rg_args_fn)
    cell_start=super_cluster_term + "_" + cell_name + "_start" +"\n"
    cell_end = super_cluster_term + "_" + cell_name + "_end" + "\n"
    string_to_write = cell_start
    grepped_transcriptome_string = subprocess.run(["rg", "-A", "1", "-f", rg_args_fn, ref_fasta_path,"|", ">>", out_fasta_fn], capture_output=True)
    grepped_transcriptome_string = grepped_transcriptome_string.stdout.decode()
    string_to_write += grepped_transcriptome_string 
    string_to_write += cell_end
    with open(out_fasta_fn, "a") as fa_out:
        fa_out.write(string_to_write)
        fa_out.close()
    subprocess.run(["echo", cell_end, ">>", out_fasta_fn]) 
    subprocess.call(["mv", out_fasta_fn, dir_name])
    subprocess.run(["mv", rg_args_fn, "last_rg_args.txt"])

