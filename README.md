# neuron_info_paper_source_code
The code for the neuron info paper

Two main goals for this paper: 

1. Create a fast method of making .fasta files containing representations of transcriptomes of all cells in an adata, that is they contain the basepair sequences for genes and their associated metadata, as well as actual lines for each read.
2. Calculate the shannon entropy of the transcriptomes created
3. Make a fast library or CLI tool to do the above, either using Cython, Rust, Mojo, C etc.
4. (if have time) Use a transformer model to create an embedding for the grammar of the total collection of reads
5. (if have time) Make a predictive model that predicts transcription changes based on exonic mutations of transcriptome sequences.


Annotated fasta files were from ensembl release 110: https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/
