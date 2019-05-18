import sys
import pandas as pd
import time
import transportation as tp

if __name__=="__main__":
    # sys.argv: 1, 销售速率文件路径； 2, 库存数据文件路径； 
    # 3, output file name;

    # sale_rate matrix, N * K, N is the number of sale nodes 
    data = pd.read_csv(sys.argv[1])
    sale_rates = data.values

    # storage data, N * K
    data = pd.read_csv(sys.argv[2])
    storages = data.values

    N = len(sale_rates)
    K = len(sale_rates[0])
    results = np.zeros((N+1,N,K))

    for k in range(K):
        sale_rate = sale_rates[:,k]
        storage = storages[:,k]
        result_k = tp.intra_alloc(sale_rate, storage)
        results[:,:,k] = result_k

    # output the result as file
    index_x = np.array([i for i in range(N+1)]) + 1
    index_y = np.array([i for i in range(N)]) + 1
    x_id = np.repeat(index_x.reshape((N+1,1)), N, axis = 1)
    y_id = np.repeat(index_y.reshape((1,N)), N+1, axis = 0)

    result_x = np.hstack((x_id.reshape(((N+1)*N,1)),y_id.reshape(((N+1)*N,1))))
    for i in range(K):
        result_x = np.hstack((result_x, results[:,:,i].reshape((N+1)*N,1)))

    data_out = pd.DataFrame(result_x)
    data_out.to_csv(sys.argv[3],index=False)
