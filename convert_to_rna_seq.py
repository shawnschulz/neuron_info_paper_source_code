#!/usr/bin/env python
import scanpy as sc 
import os
import subprocess
from queue import Queue
from threading import Thread
from time import time
import logging
import random
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-a", "--adata_name", dest="adata_name")
parser.add_option("-d", "--data_dir", dest="data_dir")
parser.add_option("-s", "--super_cluster", dest="super_cluster")
parser.add_option("-n", "--num_samples", dest="num_samples", help="supply number of cells to sample from h5ad, make sure this number is not higher than obs_names index!")

(options, args) = parser.parse_args()

data_dir = options.data_dir
adata_name = options.adata_name
super_cluster = options.super_cluster
num_samples = int(options.num_samples)

#rationale behind threading over multiprocessing is that this is an IO bound task
#want to save on memory over CPU power
#(this is actually concurrency vs parallel tasks lol)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
os.chdir(data_dir)
adata = sc.read_h5ad(adata_name)
if num_samples > len(adata.obs_names):
    raise IndexError("Please supply a number of samples that is less than or equal to the total!")

super_cluster_term = super_cluster 

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

    #testing on just 10 cells for now

obs_names = list(adata.obs_names)
random.shuffle(obs_names)
adata = adata[obs_names[0:num_samples]]

ref_fasta_path = "ensembl_annotated_transcriptome.fa"

dir_name = super_cluster_term + "_fastas"
subprocess.call(["mkdir", "-p", dir_name])


string_to_write = "" #this string is gonna be REALLY long...

cell_iteration_count=0
gene_index_chosen_count=0

#subprocess.run(["echo", "", ">", out_fasta_fn])

def append_expressed_transcripts_for_cell(cell_index):
    cell_name = adata.obs_names[cell_index]

    out_fasta_fn =  super_cluster_term + "_" + cell_name + '.fa' 
    rg_args = ""
    for gene_index in range(0, len(adata.var_names)):
        read_count = int(adata.X[cell_index, gene_index])
        ensembl_id = "gene:" + adata.var_names[gene_index]
        #should cross check the p-value of expression of the gene, hopefully it's already in the adata
        #we can use that to figure out how likely the gene is to actually be differentially expressed 
        if read_count > 0:
            rg_args = rg_args + ensembl_id + '\n'  

    rg_args_fn = cell_name + "_temp_rg_args.txt"
    with open(rg_args_fn, "w") as file:
        file.write(rg_args)
        file.close()
    cell_start=super_cluster_term + "_" + cell_name + "_start" +"\n"
    cell_end = super_cluster_term + "_" + cell_name + "_end" + "\n"
    string_to_write = cell_start
    grepped_transcriptome_string = subprocess.run(["rg", "-A", "1", "-f", rg_args_fn, ref_fasta_path], capture_output=True)
    grepped_transcriptome_string = grepped_transcriptome_string.stdout.decode()
    string_to_write += grepped_transcriptome_string 
    string_to_write += cell_end
    with open(out_fasta_fn, "a") as fa_out:
        fa_out.write(string_to_write)
        fa_out.close()
    subprocess.run(["echo", cell_end, ">>", out_fasta_fn]) 
    subprocess.call(["mv", out_fasta_fn, dir_name])
    subprocess.run(["mv", rg_args_fn, "last_rg_args.txt"])



#should spin up threads and just run all of the cell_indexes in a parallelized queue
#lets refactor to make this happen

#adata is loaded
#take the whole list of cell indexes
#run the process on multiple cells at a time, not reloading adata every time.
#for cell_index in range(0, len(adata.obs_names)): 

class transcript_worker(Thread):
    def __init__(self, queue):
        Thread.__init__(self)
        self.queue = queue
    def run(self):
        while True:
            cell_index = self.queue.get()
            try:
                append_expressed_transcripts_for_cell(cell_index)
            finally:
                self.queue.task_done()

ts = time()

queue = Queue()

#SHOVE MORE THREADS IN
for thread in range(14):
    worker = transcript_worker(queue)
    worker.daemon = True
    worker.start()

for cell_index in range(0, len(adata.obs_names)):
    queue.put(cell_index)

queue.join()
