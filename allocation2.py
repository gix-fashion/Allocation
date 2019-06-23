# -*- coding: utf-8 -*-
import sys
import pandas as pd
import time
import transportation as tp
import numpy as np

def main(argv):
    # sys.argv: 1, surplus 2, lack 3, costs 
    # 4, output file name; 5, solve method; 
    # 6, max surplus 
    
    # costs matrix, N * N, N is the number of nodes 
    # including sale nodes and warehouse nodes
    data = pd.read_csv(argv[3], index_col=0)
    costs = data.values

    # sale_pd, N * K, K is the number of cloth type
    data = pd.read_csv(argv[1], index_col=0)
    surplus = data.values[0:100]

    # storage data, N * K
    data = pd.read_csv(argv[2], index_col=0)
    lacks = data.values[0:100]

    sum_lack = np.sum(lacks, axis=1)
    print('total lack: ' + str(np.sum(sum_lack)) + ', ave_lack: ' + str(
        np.sum(sum_lack) / np.sum(np.array(sum_lack > 0).astype(int))
    ))

    sum_surplus = np.sum(surplus, axis=1)
    print('total surplus: ' + str(np.sum(sum_surplus)) + ', ave_surplus: ' + str(
        np.sum(sum_surplus) / np.sum(np.array(sum_surplus > 0).astype(int))
    ))

    if argv[5] == '0':
        result = tp.iter_trans(costs, surplus, lacks, batch = 200)
    elif argv[5] == '1':
        result = tp.intp_trans(costs, surplus, lacks, batch = 200)
    elif argv[5] == '2':
        result = tp.PSO_trans(costs, surplus, lacks, batch = 200)
    elif argv[5] == '3':
        result = tp.GA_trans(costs, surplus, lacks, batch = 200)
    '''
    N = len(surplus)
    K = len(surplus[0])

    # remove the small bags and consider the max surplus constraits
    data = pd.read_csv(argv[6], index_col=0)
    max_surp = data.values
    X_edge = np.sum(result, axis=2)
    pen_edge = np.array(X_edge < 3).astype(int)
    X_edge = X_edge * pen_edge
    result = result * np.repeat(pen_edge, K, axis=2)
    for i in range(N):
        sum_s = np.sum(X_edge[i,:])
        if sum_s > max_surp[i]:
            trans_zip = []
            for j in range(N):
                if X_edge[i][j] > 0:
                    trans_zip.append((X_edge[i][j], j))
            sort_s = sorted(trans_zip, key=lambda tzip: tzip[0])
            for sss in sort_s:
                sum_s -= sss[0]
                result[i,sss[1],:] = 0  
                if sum_s < max_surp[i]:
                    break 

    # output the result as file
    index = np.array([i for i in range(N)]) + 1
    x_id = np.repeat(index.reshape((N,1)), N, axis = 1)
    y_id = np.repeat(index.reshape((1,N)), N, axis = 0)

    result_x = np.hstack((x_id.reshape((N*N,1)),y_id.reshape((N*N,1))))
    for i in range(K):
        result_x = np.hstack((result_x, result[:,:,i].reshape(N*N,1)))

    data_out = pd.DataFrame(result_x)
    data_out.to_csv(sys.argv[4],index=False)
    '''

if __name__=="__main__":
    main(sys.argv)

