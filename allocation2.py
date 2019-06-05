# -*- coding: utf-8 -*-
import sys
import pandas as pd
import time
import transportation as tp
import numpy as np

def main(argv):
    # sys.argv: 1, surplus 2, lack 3, costs 
    # 4, output file name; 5, solve method;  
    
    # costs matrix, N * N, N is the number of nodes 
    # including sale nodes and warehouse nodes
    data = pd.read_csv(argv[3], index_col=0)
    costs = data.values

    # sale_pd, N * K, K is the number of cloth type
    data = pd.read_csv(argv[1], index_col=0)
    surplus = data.values[:,0:300]

    # storage data, N * K
    data = pd.read_csv(argv[2], index_col=0)
    lacks = data.values[:,0:300]

    if argv[5] == '0':
        result = tp.iter_trans(costs, surplus, lacks, batch = 200)
    elif argv[5] == '1':
        result = tp.intp_trans(costs, surplus, lacks, batch = 200)
    elif argv[5] == '2':
        result = tp.PSO_trans(costs, surplus, lacks, batch = 200)
    elif argv[5] == '3':
        result = tp.GA_trans(costs, surplus, lacks, batch = 200)
    
    N = len(surplus)
    K = len(surplus[0])
    # output the result as file
    index = np.array([i for i in range(N)]) + 1
    x_id = np.repeat(index.reshape((N,1)), N, axis = 1)
    y_id = np.repeat(index.reshape((1,N)), N, axis = 0)

    result_x = np.hstack((x_id.reshape((N*N,1)),y_id.reshape((N*N,1))))
    for i in range(K):
        result_x = np.hstack((result_x, result[:,:,i].reshape(N*N,1)))

    data_out = pd.DataFrame(result_x)
    data_out.to_csv(sys.argv[4],index=False)

if __name__=="__main__":
    main(sys.argv)

