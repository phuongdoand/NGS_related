#!/usr/bin/env python

import pandas as pd

mRNA_data = pd.read_csv('./mRNA_fpkm.txt', sep='\t')
mRNA_data.index = mRNA_data['Gene']
mRNA_data.drop('Gene', axis=1, inplace=True)

adj_matrix = abs(mRNA_data.T.corr('spearman'))
adj_matrix = adj_matrix.replace(1, 0)

adj_matrix.to_csv('./adj_matrix_pos.csv', index=False)
