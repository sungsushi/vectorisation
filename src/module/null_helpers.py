import numpy as np
import pandas as pd
import scipy as sp
from .fvec import bipartite_cooarray, csr_row_norm


def rn_RCshuff_wrapper(seed, a_RI_coo, meta_df, recipe_list, cuisine_list, ingredient_col):
    '''Wrapper function to recalculate the ingredient-space cuisine vectors 
    and the cuisine space ingredient vectors, from keeping Recipe-Ingredient adjacency matrix
    constant. 

    Shuffles the labels in the recipe-cuisine df meta_df.

    '''
    np.random.seed(seed)
    meta_df_bpt_shuff = meta_df.copy(True)
    type_arr = meta_df_bpt_shuff['cuisine'].to_numpy(copy=True)
    np.random.shuffle(type_arr)
    meta_df_bpt_shuff['cuisine'] = type_arr

    a_RC_coo, recipe_row, cuisine_col = bipartite_cooarray( \
        df=meta_df_bpt_shuff.sort_values(['r_id', 'cuisine']), \
        row_col=['r_id', 'cuisine'], \
        weight=False, \
        row_order=list(recipe_list), \
        col_order=list(cuisine_list))

    a_CI_csr = (a_RC_coo.T @ a_RI_coo).tocsr() # csr array

    # normalisation:
    a_CI_csr_normalised = csr_row_norm(a_CI_csr)
    norm_is_c_vec_df = pd.DataFrame(a_CI_csr_normalised.toarray(), columns=ingredient_col, index=cuisine_col)

    a_IC_csr = a_CI_csr.T
    # standarisation - correct for the number of recipes in each cuisine:
    r_sums = np.array(a_RC_coo.tocsr().sum(axis=0)).flatten()  
    r_sums[r_sums == 0] = 1 # divide by 1 instead of zero
    inv_r_sums = sp.sparse.diags(1 / r_sums)
    standardised_a_IC_csr = a_IC_csr @ inv_r_sums

    # normalisation:
    a_IC_csr_s_normalised = csr_row_norm(standardised_a_IC_csr)
    cnorm_cs_i_bpt_df = pd.DataFrame(a_IC_csr_s_normalised.toarray(), columns= cuisine_col, index=ingredient_col)
    return norm_is_c_vec_df, cnorm_cs_i_bpt_df


def rn_RIshuff_wrapper(seed, a_RC_coo, rn_df, recipe_list, cuisine_list, ingredient_list):
    '''Wrapper function to recalculate the ingredient-space cuisine vectors 
    and the cuisine space ingredient vectors, from keeping recipe-cuisine adjacency matrix
    constant. 

    Shuffles the labels in the recipe-ingredient dataframe rn_df 

    '''
    np.random.seed(seed)
    rn_df_shuff = rn_df.copy(True)
    type_arr = rn_df_shuff['r_id'].to_numpy(copy=True)
    np.random.shuffle(type_arr)
    rn_df_shuff['r_id'] = type_arr
    rn_df_shuff.drop_duplicates(inplace=True)

    a_RI_coo, recipe_row, ingredient_col = bipartite_cooarray( \
        df=rn_df_shuff.sort_values(['r_id', 'ingredient']), \
        row_col=['r_id', 'ingredient'], \
        weight=False, \
        row_order=list(recipe_list), \
        col_order=list(ingredient_list))
    # t1 = time.perf_counter()

    a_CI_csr = (a_RC_coo.T @ a_RI_coo).tocsr() # csr array

    # normalisation:
    a_CI_csr_normalised = csr_row_norm(a_CI_csr)
    norm_is_c_vec_df = pd.DataFrame(a_CI_csr_normalised.toarray(), columns=ingredient_col, index=cuisine_list)

    a_IC_csr = a_CI_csr.T
    # standarisation - correct for the number of recipes in each cuisine:
    r_sums = np.array(a_RC_coo.tocsr().sum(axis=0)).flatten()  
    r_sums[r_sums == 0] = 1 # divide by 1 instead of zero
    inv_r_sums = sp.sparse.diags(1 / r_sums)
    standardised_a_IC_csr = a_IC_csr @ inv_r_sums

    # normalisation:
    a_IC_csr_s_normalised = csr_row_norm(standardised_a_IC_csr)
    cnorm_cs_i_bpt_df = pd.DataFrame(a_IC_csr_s_normalised.toarray(), columns= cuisine_list, index=ingredient_col)
    # t2 = time.perf_counter()
    return norm_is_c_vec_df, cnorm_cs_i_bpt_df



def fw_typeshuff_wrapper(seed, meta_df, ind_to_id, a_AA_coo):
    '''Wrapper function to recalculate the animal-type vectors 
    We keep a_AA_coo constant (i.e. the food web structure) and shuffle the animal types.
    
    '''
    np.random.seed(seed)
    meta_finest_shuff = meta_df.copy(True)
    type_arr = meta_finest_shuff['type'].to_numpy(copy=True)
    np.random.shuffle(type_arr)
    meta_finest_shuff['type'] = type_arr

    shuff_bpt_AS_coo, animal_row, type_col = bipartite_cooarray(df=meta_finest_shuff, row_col=['node', 'type'], weight=False, row_order=ind_to_id)
    shuff_a_AS_out = (a_AA_coo @ shuff_bpt_AS_coo).tocsr() # out matrix
    shuff_a_AS_in = (a_AA_coo.T @ shuff_bpt_AS_coo).tocsr() # in matrix

    # out matrix normalisation:
    shuff_a_AS_out_normalised = csr_row_norm(shuff_a_AS_out)

    # in matrix normalisation:
    shuff_a_AS_in_normalised = csr_row_norm(shuff_a_AS_in)

    # put together into one dataframe:
    all_col_names = np.concatenate([type_col + '_out', type_col + '_in'])
    shuff_all_norm_vec_df = pd.DataFrame(sp.sparse.hstack([shuff_a_AS_out_normalised, shuff_a_AS_in_normalised]).toarray(), columns=all_col_names, index=animal_row)
    return shuff_all_norm_vec_df, meta_finest_shuff