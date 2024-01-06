import numpy as np
from sklearn.cluster import KMeans
import numpy.random

def convert_to_D(M, method='matrix'):
    rows = M.shape[0]
    D = np.zeros((rows, rows))
    addup = 0
    for i in range(rows):
        for j in range(i+1, rows):
            if method == 'matrix':
                D[i][j] = D[j][i] = np.sqrt(np.sum((M[i] - M[j]) ** 2))
            elif method == 'value':
                addup += np.sqrt(np.sum((M[i] - M[j]) ** 2))
    if method == 'matrix':
       return D
    elif method == 'value':
        return addup

def separate_array_by_index(arr, index_array):
    grouped_elements = {}
    for i, index in enumerate(index_array):
        if index not in grouped_elements:
            grouped_elements[index] = []
        grouped_elements[index].append(arr[i])

    # Convert the dictionary values to a list of lists
    separated_arrays = list(grouped_elements.values())
    return separated_arrays


def find_k(M):
    gap_list = []
    for k in range(2, 8):
        D = []
        W = []
        kmeans = KMeans(n_clusters = k, random_state=0).fit(M)
        index_array = kmeans.labels_
        separated_arrays = separate_array_by_index(M_copy, index_array)
        for r in range(len(separated_arrays)):
            D[r] = convert_to_D(separated_arrays[r], method='value')
            W[r] = D[r] / (2 * len(separated_arrays[r]))
        W_sum = np.sum(W)
        W_sum_log = np.log(W_sum)

        # Reference distribution
        ref_Wk = []
        for b in range(10):
            # Create reference dataset
            M_ref = np.random.random_sample(M.shape)
            kmeans = KMeans(n_clusters=k, random_state=b).fit(M_ref)
            index_array = kmeans.labels_
            separated_arrays = separate_array_by_index(M_ref, index_array)
            Wk_ref = np.sum([convert_to_D(cluster, method='value') / (2 * len(cluster)) 
                             for cluster in separated_arrays])
            ref_Wk.append(np.log(Wk_ref))
        W_ref_log = np.mean(ref_Wk)

        # Calculate the gap statistic
        gap = W_ref_log - W_sum_log
        gap_list.append(gap)

    best_k = 0
    for i in range(1, len(gap_list)):
        if gap_list[i] > gap_list[k]:
            best_k = i
    print(best_k + 2)   



#Test code
matrix = np.array([
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9],
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9],
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9]
])

find_k(matrix)