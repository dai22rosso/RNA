import numpy as np
import sympy as sp
import anndata as ad
import pandas as pd
adata=ad.read_h5ad('./pancreas.h5ad')
adata.X.todense()#将稀疏矩阵转成普通矩阵
X=pd.DataFrame(adata.X.todense())
cell_name=adata.obs.index
print('cell',cell_name,len(cell_name))
chr_name=adata.var.index
print('chr',chr_name,len(chr_name))
X.index=cell_name
X.columns=chr_name
X=np.array(X.T) #行为peak,列为cell
print(X[0,0],'X',X.shape)
#X.to_csv('./RNA_pancreas.csv',sep='\t')
