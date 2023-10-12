#!/usr/bin/env python
# coding: utf-8

# # Analyze the BD results

# In[1]:


import scvelo as scv
import scanpy
import igraph
import glob, os
import pandas as pd
import numpy as np
import re
from collections import Counter
import anndata

import shutil
    
import h5py
from shutil import copyfile


def copyFiles(files, to):
    for f in files:
        name = os.path.basename( f )
        print( f"copy {f} to {to}" )
        copyfile( f, os.path.join(to, name ) )
    print( "all copied" )


# ## This now (2023.05.20)
# 
# is based on the Rhapsody 1.2.1 (first commit).
# The data was created in the TestData folder using the commands
# 
# ```
# ../target/release/quantify_rhapsody -r cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz -f cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz -o BD_results/Rustody_S1 -s mouse  -e 2276_20220531_chang_to_rpl36a_amplicons.fasta -a MyAbSeqP
# anel.fasta -m 200 -v v2.96 --gene-kmers 16
# 
# ../target/release/quantify_rhapsody -r cells.1.Rhapsody_SV_index2_S2_R1_001.fastq.gz -f cells.1.Rhapsody_SV_index2_S2_R2_001.fastq.gz -o BD_results/Rustody_S2 -s mouse  -e 2276_20220531_chang_to_rpl36a_amplicons.fasta -a MyAbSeqPanel.fasta -m 200 -v v2.96 --gene-kmers 16
# ```
# 
# The --gene_kmers makes a huge difference - still.
# This is the setting where I hat more antibody tags out.
# 
# I am shure this mapper is more stringent than my first implementation as this here actually checks two consecutive DNA fragments 8bp +16bp. Hence the matching are is bigger than in version 0.x. 
# 
# ## Overall result
# 
# Looks acceptable. I identify less reads as containing sequence, but the result still looks more or less the same.
# I need to look into that in more detail.

# In[2]:


print('\n'.join(f'{m.__name__}=={m.__version__}' for m in globals().values() if getattr(m, '__version__', None)))


# In[3]:


files = [ 'Combined_S1Subset_DBEC_MolsPerCell.csv', 'Combined_S2Subset_DBEC_MolsPerCell.csv' ]
sinfoF = [ 'S1Subset_Sample_Tag_Calls.csv', 'S2Subset_Sample_Tag_Calls.csv']
sname = ['S1_BD', 'S2_BD' ]

ofile = "BD_compared.h5ad"
ofile_raw = "BD_compared_raw.h5ad"

# In[4]:


df = pd.read_csv(files[0], skiprows=7, index_col=0)
df.columns[0:20]


# In[5]:


m = np.array([ s.__contains__ ("|") for s in df.columns ]).sum()


# In[6]:


adata = anndata.AnnData(X = df.iloc[0:,m:])


# In[7]:


adata.var_names
adata.obs['sample'] = sname[0]


# In[8]:


df = pd.read_csv(files[1], skiprows=7, index_col=0)
tmp = anndata.AnnData(X = df.iloc[0:,m:])
tmp.obs['sample'] = sname[1]
adata = adata.concatenate( tmp )


# In[9]:


os.system('  cargo build -r  & source ~/.cargo/env')


# In[10]:


os.name


# In[18]:


if os.path.exists("Rustody_S1"):
    shutil.rmtree('Rustody_S1')
if not os.path.exists("Rustody_S1"):
    f1 =  os.path.join( "..", "cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz")
    f2 =  os.path.join( "..", "cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz")
    ab =  os.path.join( "..", "MyAbSeqPanel.fasta")
    ex =  os.path.join( "..", "2276_20220531_chang_to_rpl36a_amplicons.fasta")
    exe = os.path.join( "..", '..', 'target', 'release', 'quantify_rhapsody_multi' )
    if os.name == 'nt': # windows...
        exe = exe + '.exe'
    print ( f"{exe} -r {f1} -f {f2} -o Rustody_S1 -s mouse  -e {ex} -a {ab} -m 200 -v v2.96")
    os.system(f"{exe} -r {f1} -f {f2} -o Rustody_S1 -s mouse  -e {ex} -a {ab} -m 200 -v v2.96")



# In[20]:


if os.path.exists("Rustody_S2"):
    shutil.rmtree('Rustody_S2')
if not os.path.exists("Rustody_S2"):
    f1 =  os.path.join( "..", "cells.1.Rhapsody_SV_index2_S2_R1_001.fastq.gz")
    f2 =  os.path.join( "..", "cells.1.Rhapsody_SV_index2_S2_R2_001.fastq.gz")
    ab =  os.path.join( "..", "MyAbSeqPanel.fasta")
    ex =  os.path.join( "..", "2276_20220531_chang_to_rpl36a_amplicons.fasta")
    print ( f"{exe} -r {f1} -f {f2} -o Rustody_S2 -s mouse  -e {ex} -a {ab} -m 200 -v v2.96")
    os.system(f"{exe} -r {f1} -f {f2} -o Rustody_S2 -s mouse  -e {ex} -a {ab} -m 200 -v v2.96")


# In[21]:


def readRustodyExpression(path, name):
    print(f"reading Rustody expression from path {path}/BD_Rhapsody_expression/")
    this = scanpy.read_10x_mtx( path+'/BD_Rhapsody_expression/' )
    this.obs['sample'] = name
    obs1 = pd.read_csv( path+'/SampleCounts.tsv', sep="\t")
    this.obs = this.obs.merge( obs1, left_index= True, right_on = 'CellID' )
    this.obs_names = this.obs['CellID'] + "_" +  this.obs['sample']
    this = this[this.obs['AsignedSampleName'] != "na"]
    this.obs['AsignedSampleName'] = this.obs['AsignedSampleName'] + "_" + this.obs['sample']
    return(this)


# In[22]:


adata = adata.concatenate( readRustodyExpression( 'Rustody_S1', 'S1_Rustody'))
adata = adata.concatenate(readRustodyExpression( 'Rustody_S2', 'S2_Rustody'))
adata.write( ofile_raw )
print(f"I have written the raw data to {ofile_raw}")

# In[23]:


Counter( adata.obs['sample'])


# In[24]:


adata.var_names.values


# In[25]:


adata[:,"Cnot2"][adata.obs['sample'] == 'S1_Rustody'].X.todense().sum()


# In[26]:


scanpy.pp.filter_genes(adata, min_counts=1 )
scanpy.pp.filter_cells(adata, min_genes=20 )
Counter( adata.obs['sample'])


# In[27]:


adata.obs['cellID'] =  [ re.findall(r'\d+', n)[0] for n in adata.obs_names] 


# In[28]:


Counter(adata.obs['cellID'])


# In[29]:


#dat = adata[ adata.obs['cellID'] == '4753157' ].X.todense()
#ids = [id  for id in range(dat.shape[1]) if not dat[0,id] == dat[1,id] ]
#bools = [ not dat[0,id] == dat[1,id]  for id in range(dat.shape[1]) ]
#adata.var.loc[bools,:]


# In[30]:


adata


# In[31]:


#diff_val_BD = [ int(n) for n in dat[0,bools].transpose()]


# In[32]:


#diff_val_Rhapsody = [ int(n) for n in dat[1,bools].transpose()]


# In[33]:


#import plotly.express as px
#df = pd.DataFrame( { 'BD' :np.log1p(diff_val_BD), 'Rhapsody' : np.log1p(diff_val_Rhapsody)})
#fig = px.violin(df, y= df.columns )
#fig.show()
#

# In[34]:


#adata[ adata.obs['cellID'] == '4753157' ].obs


# In[35]:


#dat.shape[1]


# In[36]:


scanpy.pp.calculate_qc_metrics( adata, inplace=True, log1p=True, percent_top= [10])


# In[37]:


scanpy.pp.filter_cells(adata, min_counts=400 )
scanpy.pp.downsample_counts(adata, counts_per_cell= 400 )
scv.pp.log1p(adata)
Counter( adata.obs['sample'])


# In[38]:


scanpy.pp.neighbors(adata)
dimensions = 2
scanpy.tl.umap(adata,n_components= dimensions)


# In[39]:


scv.pl.scatter(adata, color='sample', figsize =(7,5), legend_loc='right margin')


# In[40]:


scanpy.tl.louvain(adata)


# In[41]:


scv.pl.scatter(adata, color='louvain', figsize =(7,5), legend_loc='right margin')


# In[42]:


adata


# In[43]:


adata.obs.pivot_table(values = "n_counts", index = "louvain", columns="sample", aggfunc='count')


# In[44]:


adata.obs['n_genes' ] = adata.obs['n_genes' ].astype('int32')


# In[45]:


adata


# In[46]:


scanpy.pl.violin(adata, [ 'n_genes', 'log1p_total_counts' ], 'sample')


# In[47]:


key_added = "louvain"
scanpy.tl.rank_genes_groups(
    adata, 
    groupby   = 'louvain',
    key_added = key_added,
    method    = 'wilcoxon',
)

scanpy.pl.rank_genes_groups(adata, key = key_added )


# In[48]:


scv.pl.scatter(adata, color=['Igkc', 'Igha', 'Ighm'], figsize =(7,5), legend_loc='right margin')


# In[49]:


scv.pl.scatter(adata, color=[ 'Ighg1', 'Ighg2b', 'Ighg3' ], figsize =(7,5), legend_loc='right margin')


# In[50]:


names = np.array(adata.obs['sample'].unique()).astype('str')
names.sort()
names


# In[51]:


for n in names:
    print (n)
    scv.pl.scatter(adata[adata.obs['sample'] == n], color=[ 'Ighg1', 'Ighg2b', 'Ighg3' ], figsize =(7,5), legend_loc='right margin')


# In[52]:


for n in names:
    print (n)
    scv.pl.scatter(adata[adata.obs['sample'] == n], color=[ 'Top2a', 'Tmem97', 'Jchain' ], figsize =(7,5), legend_loc='right margin')


# In[ ]:


adata.obs['CellID'] = adata.obs_names
adata.obs['AsignedSampleName'] = 'unimportant'
adata.write(ofile)
print ( f"data written to {ofile}")

# In[ ]:





# In[ ]:




