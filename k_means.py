import numpy as np
import sympy as sp
import anndata as ad
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix

def get_index1(lst=None, item=''):
    tmp = []
    tag = 0
    for i in lst:
        if i == item:
            tmp.append(tag)
        tag += 1
    return tmp

def coo_matrix_transform(A):
    A_sparse = coo_matrix(A)
    output = []
    for i in range(0, len(A_sparse.data)):
        output.append([A_sparse.row[i], A_sparse.col[i], A_sparse.data[i]])
    return output
def find_dii1(M1, M2):
    # Data Structure of M1/M2: arrays of [col, value]
    # M1/M2 are sorted by col
    M1 = sorted(M1, key=lambda x: x[0])
    M2 = sorted(M2, key=lambda x: x[0])
    pointer1 = 0
    pointer2 = 0
    diff = 0
    while pointer1 < len(M1) or pointer2 < len(M2):
        if pointer1 == len(M1):
            diff += M2[pointer2][1] ** 2
            pointer2 += 1
        elif pointer2 == len(M2):
            diff += M1[pointer1][1] ** 2
            pointer1 += 1
        elif M1[pointer1][0] == M2[pointer2][0]:
            diff += (M1[pointer1][1] - M2[pointer2][1]) **2
            pointer1 += 1
            pointer2 += 1
        elif M1[pointer1][0] < M2[pointer2][0]:
            diff += M1[pointer1][1] ** 2
            pointer1 += 1
        else:
            diff += M2[pointer2][1] ** 2
            pointer2 += 1
    return diff
def matrix_dk(A,k):
    A2=coo_matrix_transform(A)
    
    d=np.zeros([len(A),len(A)])
    for i in range(len(A)):
        li=get_index1(A2[:][0],i)
        for i1 in range(len(A)):
            if(i>i1):
                continue
            else:
                
                li1=get_index1(A2[:][0],i1)
                M1=[]
                M2=[]
                for k in range(len(li)):
                    M1.append(A2[li[k]])
                for k in range(len(li1)):
                    M2.append(A2[li1[k]])
                dii1=find_dii1(M1, M2)
                d[i][i1]=dii1
                d[i1][i]=dii1
    return d

def d_to_D(d,A,k):
    clustering = KMeans(n_clusters=k, random_state=1).fit_predict(A)
    D=np.zeros(k)
    lD=np.zeros(k)
    for i in range(k):
        cri=get_index1(clustering,i)
        sum_d=0
        for j in range(len(cri)):
            for j1 in range(len(cri)):
                sum_d+=d[j][j1]
        lD[i]=len(cri)
        D[i]=sum_d
    return D,lD
                    
                    

def D_to_Wk(D,A, k,lD):
    # Ensure A is a sparse CSR matrix
    Wk=0
    for i in range(k):
        Wk+=D[i]/lD[i]
   

    return Wk

def matrix_to_Wk(A,k):
    d=matrix_dk(A, k)
    D,lD=d_to_D(d, A, k)
    return D_to_Wk(D, A, k, lD)

#X:matrix to be dealed with
#n:the most 1/n gene will be picked
def pick_gene_piece(X,n):
    X_1=np.zeros((int(len(X)/n),len(X[0])))
    sum=np.zeros(len(X))
    for i in range(len(X)):
        sum[i]=np.sum(X[i])
    rank=list(np.argsort(sum))
    for i in range(int(len(X)/n)):
        X_1[i]=X[rank.index(i)]
    return X_1

def sort_cell(sort_num,pred,cell_name):
    A=[]
    for i in range(sort_num):
        l=[]
        for ii in range(len(pred)):
            if(pred[ii]==i):
                l.append(cell_name[ii])
            else:
                continue
        A.append(l)
    return A

def K_estimate(X,k):
    max1=0
    for i in range(len(X)):
        if(max(X[i])>max1):
            max1=max(X[i])
        else:
            continue
    print('pass max',k)
    lg=np.zeros(20)
    for i in range(20):
        X1=X.copy()
        for ii in range(len(X)):
            n=np.random.random(len(X[0]))
            n=n*max1
            X1[ii]=n
        print('kazhule',k)
        Wk=matrix_to_Wk(X1,k)
        lg[i]=np.log10(Wk)
    print('pass matrix',k)
    Wkmean=np.mean(lg)
    G=Wkmean-matrix_to_Wk(X,k)
    return G

def get_k(X):
    G=np.zeros(30)
    for k in range(30):
        print('enter loop',k)
        G[k]=K_estimate(X,k+1)
        
    G=list(G)
    return G.index(max(G))+1

    
    
        
        
        
    
    

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
X2=pick_gene_piece(X,100)
X2=X2.T
pca=PCA(n_components=30)
pca.fit(X2)
M=pca.transform(X2)
print(M.shape)
x=np.zeros(len(M))
y=x.copy()
z=x.copy()
for i in range(len(M)):
    y[i]=M[i][1]
    x[i]=M[i][0]
    # z[i]=M[i][2]
print(x,y)
plt.scatter(x,y,s=0.1)
plt.show()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# # Set the size of the dots (adjust the 's' parameter)
# dot_size = 30

# ax.scatter(x, y, z, s=0.008)

# Add labels and title
# ax.set_xlabel('X-Axis')
# ax.set_ylabel('Y-Axis')
# ax.set_zlabel('Z-Axis')
# plt.title('3D Scatter Plot')

# Show the plot
#plt.show()
k=get_k(M)
y_pred = KMeans(n_clusters=k, random_state=1).fit_predict(M)
A=sort_cell(k,y_pred,cell_name)
A=np.array(A)
print(A.shape)
