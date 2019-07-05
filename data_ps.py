import pandas as pd
import numpy as np
import math
import sys
import os

supply = pd.read_csv('./data/supply_0623.csv', encoding='gbk')

pre_results = pd.read_csv('./data/623_results_nm.csv', encoding='utf-8')
data_len = len(pre_results.values)

data = pd.read_csv('./data/trans_time.csv')
provinces = data.values[:,0]
ptimes = data.values[:,1:]
# generate province dict
pro_id = [i for i in range(len(provinces))]
dic_pro = dict(zip(provinces,pro_id))

# generate store dict
list_store = list(set(pre_results['ShopCode']))
N = len(list_store)
store_id = [i for i in range(N)]
dic_store = dict(zip(list_store, store_id))

# generate skc dict
SKC_data = pre_results['MatCode'] + pre_results['SizeName']
list_skc = list(set(SKC_data))
K = len(list_skc)
skc_id = [i for i in range(K)]
dic_skc = dict(zip(list_skc, skc_id))
skc_mat_size = []
for i in range(K):
    sku_pos = np.where(SKC_data == list_skc[i])[0][0]
    skc_mat_size.append([pre_results['MatCode'][sku_pos],pre_results['SizeName'][sku_pos]])

# generate map from store to province
store2pro = np.zeros((N,))
for i in range(N):
    s_poss = np.where(supply['ShopCode'] == list_store[i])
    if len(s_poss[0]) == 0:
        print(list_store[i])
        s_pos = np.random.randint(len(supply['ShopCode']))
    else:
        s_pos = s_poss[0][0]
    s2p = supply['ProvinceName'][s_pos][0:2]
    s_pro = dic_pro[s2p]
    if s_pro < 0 or s_pro >= len(provinces):
        print(supply['ProvinceName'][i][0:2])
    store2pro[i] = s_pro

# generate trans time between stores
store2pro = store2pro.astype(int)
stimes = np.zeros((N,N))
for i in range(N):
    stimes[i,:] = ptimes[store2pro[i]*np.ones((1,N)).astype(int),store2pro]

print('store num: ' + str(N) + ', skc num: ' + str(K))

# handle needs
surplus = np.zeros((N,K))
lack = np.zeros((N,K))
for i in range(data_len):
    if pre_results['SkuNeed'][i] > 0:
        lack[dic_store[pre_results['ShopCode'][i]]][dic_skc[SKC_data[i]]] = pre_results['SkuNeed'][i]
    elif pre_results['SkuNeed'][i] < 0:
        surplus[dic_store[pre_results['ShopCode'][i]]][dic_skc[SKC_data[i]]] = - pre_results['SkuNeed'][i]

store_Frame = pd.DataFrame(np.array(list_store))
store_Frame.to_csv('~/production/data_623/store_list_nm.csv') # the store list, 其顺序决定了allocation2.py输出时店铺的编号

sku_Frame = pd.DataFrame(skc_mat_size)
sku_Frame.to_csv('~/production/data_623/sku_list_nm.csv')   # sku list, 同上

surp_Frame = pd.DataFrame(surplus)
surp_Frame.to_csv('~/production/data_623/surplus_nm.csv')   # 供给量csv文件

lack_Frame = pd.DataFrame(lack)
lack_Frame.to_csv('~/production/data_623/lack_nm.csv')      # 需求量csv文件

cost_Frame = pd.DataFrame(stimes)
cost_Frame.to_csv('~/production/data_623/costs_nm.csv')     # 运输时间csv文件

print(1)
