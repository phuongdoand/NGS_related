#!/usr/bin/env python

import pandas as pd
import numpy as np
import random
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-a", help = "Adjacency matrix file", dest = "a")
parser.add_argument("-s", help = "Seed genes list", dest = "s")
parser.add_argument("-o", help = "Result of permutation test", dest = "o")
args = parser.parse_args()

def power_fun(k, r, e):   
    return r*k + e

adj_matrix = pd.read_csv(args.a)

# Create the weighted co-expression network
A_mat = np.matrix((adj_matrix**4).values)
for i in range(len(A_mat)):
    A_mat[i, :] = A_mat[i, :] / np.sum(A_mat[i, :])
W_mat = A_mat.T

# RWR funcion
r = 0.7

def RWR(r, P0, W):
    i = 0
    Pi = P0
    while i >= 0:
        Pi_1 = (1 - r)*W*Pi + r*P0
        if np.linalg.norm(Pi_1 - Pi) < 10**(-6):
            break
        else:
            Pi = Pi_1
            i += 1
    return Pi_1

# Read the seed genes and prepare for RWR
seed_genes = pd.read_csv(args.s, sep='\t')
seed_genes = pd.DataFrame(seed_genes['Genes'])
genes_df = pd.DataFrame({'Genes': adj_matrix.columns.values})
genes_df.reset_index(inplace=True)
seed_index = pd.merge(genes_df, seed_genes, on='Genes')

# construct P0 matrix (put 1/n into P0 through the index of seed genes)
P0_seed = np.zeros((len(genes_df), 1))
P0_seed[seed_index.loc[:, 'index']] = 1 / len(seed_index)
P0_seed = np.matrix(P0_seed)

# RWR
Pi_mat = RWR(r, P0_seed, W_mat)
result = pd.concat([genes_df, pd.DataFrame({'probability': np.array(Pi_mat).flatten().tolist()})], axis=1)
result.drop(seed_index['index'], inplace=True)
rank = result.drop('index', axis=1)
rank = rank.sort_values('probability', ascending=False).reset_index(drop=True).reset_index()
rank['index'] = rank['index'] + 1
rank.rename(columns={'index': 'ranking'}, inplace=True)

# Permutation test
# construct the real RWR result for comparison
true = pd.DataFrame({'probability': np.array(Pi_mat).flatten().tolist()})

# construct random P0 matrix, perform permutation test
i = 0
count_no = np.zeros((len(genes_df), 1))
while i < 1000:
    random.seed(a=i)
    random_index = random.sample(list(genes_df.loc[:, 'index']), len(seed_index))
    P0_random = np.zeros((len(genes_df), 1))
    P0_random[random_index] = 1 / len(random_index)
    P0_random = np.matrix(P0_random)

    Pi_random = pd.DataFrame({'probability': np.array(RWR(r, P0_random, W_mat)).flatten().tolist()})
    Pi_random.iloc[random_index,:] = -1  # avoid genes that being as seeds, larger than the real probability

    r_vs_t = (Pi_random > true)['probability']
    count_no[list(r_vs_t[r_vs_t == True].index)] = count_no[list(r_vs_t[r_vs_t == True].index)] + 1
    i += 1


genes_pvalue = pd.concat([genes_df, pd.DataFrame({'P-value': (count_no / 1000).flatten().tolist()})], axis=1)
genes_pvalue.drop(seed_index['index'], inplace=True)
genes_pvalue.drop('index', axis=1, inplace=True)

rank_plus = pd.merge(rank, genes_pvalue, on='Genes')
rank_plus_5out = rank_plus[rank_plus['P-value'] < 0.05]
rank_plus_5out.reset_index(drop=True, inplace=True)

# Save the result of permutation test
rank_plus_5out.to_csv(args.o, index=False)



