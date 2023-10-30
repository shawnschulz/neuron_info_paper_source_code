#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
from optparse import OptionParser


# In[ ]:


def prepare_options(parser):
    """Prepare options parser
    """
    parser.add_option("-i", "--input", dest="input_file",
                      help="Input FILE", metavar="FILE")
    parser.add_option("-n", "--name", dest="name",
                      help="filename", metvar="FILENAME")
(options, args) = parser.parse_args()
adata_path = options.input_file
adata_name = options.name


# In[ ]:


ts_adata = sc.read_h5ad(adata_path)


# In[ ]:


sc_ts_all_qc_df = sc.pp.calculate_qc_metrics(ts_adata)


# In[ ]:


var_qc_df =sc_ts_all_qc_df[1]


# In[ ]:


obs_qc_df = sc_ts_all_qc_df[0]


# In[ ]:


obs_qc_df["supercluster_term"] = ts_adata.obs["supercluster_term"]


# In[ ]:


obs_qc_df.to_csv("../" + adata_name + "_obs_qc.csv")


# In[ ]:


var_qc_df["gene_id"] = ts_adata.var['feature_name']


# In[ ]:


var_qc_df.to_csv("../" + adata_name + "ts_all_var_qc.csv")


# In[ ]:


n_genes_df = ts_adata.obs[['total_genes', 'supercluster_term']]


# In[ ]:


total_genes_by_cell_type = n_genes_df.groupby('cell_type').mean().T.sort_values('total_genes', axis=1, ascending=False)


# In[ ]:


total_genes_by_cell_type.to_csv("../" + adata_name + "ts_all_tot_genes_by_cell_type.csv")

