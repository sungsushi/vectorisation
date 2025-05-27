import numpy as np
from scipy.cluster.hierarchy import to_tree
import pandas as pd
import scipy as sp 
import matplotlib.pyplot as plt

def get_clustersize_range(Z, ind_to_id, sizes):
    '''Get the clusters of at least size members. '''
    
    rootnode, nodelist = to_tree(Z, rd=True)
    bools = [i.is_leaf() for i in nodelist]
    opp_bools = list(map(lambda x: not x, bools))
    # leaves_are = np.array([i.get_id() for i in nodelist])[bools]
    # not_leaves = np.array([i.get_id() for i in nodelist])[opp_bools]

    # preorder == implicit that the clusters are all leaves! They're consistuent leaves...
    lst = [i.pre_order(lambda x: ind_to_id[x.id]) for i in np.array(nodelist)[opp_bools]] # getting list of all clusters 
    clstr_sizes = [len(i) for i in lst] # sizes of clusters
    cluster_bool = (np.array(clstr_sizes) >= sizes[0])*(np.array(clstr_sizes) <= sizes[1] )  # only consider clusters of at least size
    dists = [i.dist for i in np.array(nodelist)[opp_bools]] # getting distances of all clusters 

    rel_clusters = np.array(lst, dtype=object)[cluster_bool] # all clusters containing clstrsz number. 
    rel_dists = np.array(dists, dtype=object)[cluster_bool] # all cluster distances containing clstrsz number.|
    rel_cluster_tuples = [tuple(i) for i in rel_clusters]

    rel_df = pd.DataFrame({'cluster':rel_cluster_tuples, f'dists':rel_dists}).sort_values(by='dists', ascending=True)
    return rel_df

def hypergeometric_pmf(k, K, n, N):
    '''
    k : number of successes in sample
    K : number of total labelled successes in population
    n : number in sample
    N : total population size

    '''

    t1 = sp.special.comb(K,k, repetition=False)
    t2 = sp.special.comb(N-K,n-k, repetition=False)
    t3 = sp.special.comb(N, n, repetition=False)

    p = (t1 * t2 / t3)

    return p


def one_test(k, n, K, N):
    '''one-tailed hypothesis test'''
    ks = range(k, n+1) # assumes extremity for larger k values
    ps = 0.

    for kay in ks:
        p_val = hypergeometric_pmf(k=kay, K=K, n=n, N=N)
        ps += p_val
        # print(f'k = {kay},', 'p =', hypergeometric_pmf(k=kay, K=K, n=n, N=N))
        
    return ps




def row_count_subset(csrarray, r_ids, c_ids, row_names, col_names):
    """
    Extracts a subset of rows and columns from a sparse matrix and computes the frequency of column IDs.

    Parameters:
    csrarray (scipy.sparse.csr_matrix): The sparse matrix to slice.
    r_ids (list): List of row IDs to include in the subset.
    c_ids (list): List of column IDs to include in the subset.
    row_names (list): List of all row names corresponding to the rows in the sparse matrix.
    col_names (list): List of all column names corresponding to the columns in the sparse matrix.

    Returns:
    dict: A dictionary where keys are column IDs (from `c_ids`) and values are their summed frequencies 
          across the selected rows (from `r_ids`). Only columns with non-zero frequencies are included.
    """
    cids_bool = pd.Series(col_names).isin(c_ids).values  # Boolean mask for columns of interest
    rids_bool = pd.Series(row_names).isin(r_ids).values  # Boolean mask for rows of interest
    rids_go_array = csrarray[rids_bool]  # Slice the sparse matrix by rows
    rids_bool_summed = rids_go_array.sum(axis=0)  # Sum over the selected rows
    sliced_go_counts = rids_bool_summed[cids_bool]  # Filter counts for the selected columns
    sliced_col_names = col_names[cids_bool]  # Get column names for the selected columns

    # Create a dictionary of column IDs and their frequencies, excluding zero-frequency entries
    ids_counter = dict(zip(sliced_col_names[sliced_go_counts > 0], sliced_go_counts[sliced_go_counts > 0]))
    return ids_counter




def make_GO_table(data, ax, highlight_row=None):

    ax.axis('off')

    # Create the table
    table = ax.table(
        cellText=data,
        colLabels=['GO terms', 'Description'],
        cellLoc='center',
        loc='center'
    )

    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1, 1.8)

    col_widths = [0.2,0.8]  # fractions of figure width
    # nrows, ncols = df.shape

    for i in range(len(data)+1):
        for j in range(2):
            cell = table[i, j]

            cell.set_width(col_widths[j])
            cell.set_edgecolor('black')


    GO_col = 0  # "Dist"
    desc_col = 1

    for row in range(len(data)+1):
        GO_col_cell = table[(row, GO_col)]  # +1 to skip header row
        text = GO_col_cell.get_text()
        text.set_fontweight('bold')
        GO_col_cell = table[(row, desc_col)]  # +1 to skip header row
        GO_col_cell.get_text().set_ha('left')
        desc_col_cell = table[(row, desc_col)]
        text = desc_col_cell.get_text().get_text()
        if '\n' in text:
            highlight_row = row

    if highlight_row is not None:
        # print(highlight_row)
        for col in range(2):  # Number of columns
            table[(highlight_row, col)].set_height(0.1) 

    # Optional: bold the column header too
    n_cols = 2
    for col in range(n_cols):
        header_cell = table[(0, col)]
        header_cell.get_text().set_fontweight('bold')    # plt.show()
        header_cell.get_text().set_fontsize(15)

def insert_line_break(s):
    if len(s) > 50:
        # Find the last space before or at the 55th character
        break_pos = s.rfind(' ', 0, 50)
        if break_pos == -1:
            # No space found before 55 — fall back to hard break
            break_pos = 50
        return s[:break_pos] + '\n' + s[break_pos+1:]
    return s
