import scanpy as sc
import scipy.io
import pandas as pd
import pickle
import sys

if __name__ == '__main__':
    input = sys.argv[1]
    output = sys.argv[2]
    with open(input, 'rb') as f:
        a=pickle.load(f)
    a.var.to_csv(output+'_var.csv')
    a.obs.to_csv(output+'_obs.csv')
    # write_matrix_data()
    