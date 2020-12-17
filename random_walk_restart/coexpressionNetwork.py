#!/usr/bin/env python

import pandas as pd

adj_matrix = pd.read_csv('adj_matrix_pos.csv')
A_mat = np.matrix((adj_matrix**6).values)



