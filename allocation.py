import sys
import pandas as pd
import time
import transportation as tp

if __name__=="__main__":
    # sys.argv: 1, 销售预测数据文件路径； 2, 库存数据文件路径； 
    # 3, 期货与补货数据文件路径； 4, 运费数据文件; 5, output 
    # file name; 6, solve method;  
    
    # costs matrix, N * N, N is the number of nodes 
    # including sale nodes and warehouse nodes
    data = pd.read_csv(sys.argv[4])
    costs_edge = data.values

    # sale_pd, N * K, K is the number of cloth type
    data = pd.read_csv(sys.argv[1])
    sale_pd = data.values

    # storage data, N * K
    data = pd.read_csv(sys.argv[2])
    storage = data.values

    # supplement data, N * K
    data = pd.read_csv(sys.argv[3])
    supp = data.values

    # compute the surplus and lacks
    total = storage + supp - sale_pd
    surplus = total * np.array(total > 0).astype(int)
    lacks = - total * np.array(total < 0).astype(int)

    if sys.argv[6] == 0:
        result = tp.iter_trans(costs, surplus, lacks)
    elif sys.argv[6] == 1:
        result = tp.intp_trans(costs, surplus, lacks)
    
    N = len(sale_pd)
    K = len(sale_pd[0])
    # output the result as file
    index = np.array([i for i in range(N)]) + 1
    x_id = np.repeat(index.reshape((N,1)), N, axis = 1)
    y_id = np.repeat(index.reshape((1,N)), N, axis = 0)

    result_x = np.hstack((x_id.reshape((N*N,1)),y_id.reshape((N*N,1))))
    for i in range(K):
        result_x = np.hstack((result_x, result[:,:,i].reshape(N*N,1)))

    data_out = pd.DataFrame(result_x)
    data_out.to_csv(sys.argv[5],index=False)

