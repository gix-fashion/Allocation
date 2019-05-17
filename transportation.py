from pulp import *
import numpy as np
import math

def trans(costs, x_max, y_max):

    row = len(costs)
    col = len(costs[0])
    Rows = range(row)
    Cols = range(col)

    prob = LpProblem("Transportation Problem", LpMinimize)
    vars = LpVariable.dicts("x",(Rows,Cols),lowBound = 0)

    # prob += lpSum([lpSum([vars[i][j] * costs[i][j] for j in Cols]) for i in Rows])
    prob += lpSum([lpDot([vars[i][j] for j in Cols], costs[i]) for i in Rows])

    for i in range(row):
        prob += lpSum([vars[i][j] for j in Cols]) == x_max[i]

    for j in range(col):
        prob += lpSum([vars[i][j] for i in Rows]) == y_max[j]

    prob.solve()
    if LpStatus[prob.status] != 'Optimal':
        print("Status:", LpStatus[prob.status])
    return {'objective':value(prob.objective), 'var': [[value(vars[i][j]) for j in range(col)] for i in range(row)]}

def muti_trans(costs, surp, lack):
    N = len(costs) # number of total area
    K = len(surp[0]) # type number of total cloth

    X = np.zeros((N,N,K))
    obj = 0
    for i in range(K):
        # target nodes
        l_k = lack[:,i]
    
        # source nodes
        s_k = surp[:,i]

        T_s = np.sum(s_k) - np.sum(l_k)

        index_l = np.where(l_k > 0)[0]
        index_s = np.where(s_k > 0)[0]

        # costs of k type
        costs_k = costs[index_s,:]
        costs_k = costs_k[:,index_l]
        costs_k = np.reshape(costs_k,(len(index_s),len(index_l)))
        # lack of k type
        lack_k = l_k[index_l]
        surplus_K = s_k[index_s]

        if T_s == 0:
            # solve transportation problem
            res = trans(costs_k, surplus_K, lack_k)
            x_k = np.array(res['var'])
            obj_k = res['objective']
        elif T_s < 0: # the surplus is less than lack
            cost_add = np.zeros((1,len(index_l)))
            costs_k = np.vstack((costs_k,cost_add))
            surplus_K = np.hstack((surplus_K,-T_s))
            res = trans(costs_k, surplus_K, lack_k)
            x_k = np.array(res['var'])
            x_k = x_k[0:len(index_s),:]
            obj_k = res['objective']
        else:
            cost_add = np.zeros((len(index_s),1))
            costs_k = np.hstack((costs_k,cost_add))
            lack_k = np.hstack((lack_k,T_s))
            res = trans(costs_k, surplus_K, lack_k)
            x_k = np.array(res['var'])
            x_k = x_k[:,0:len(index_l)]
            obj_k = res['objective']

        for li in range(len(index_l)):
            X[index_s,index_l[li],i] = x_k[:,li]
        
        obj += obj_k

    return {'objective':obj, 'var': X}

def updata_costs(costs, batch, x_edge, stair, p1 = 5, p2 = 1):
    if stair == 0:
        bags = np.ceil(x_edge / batch)
        margin = x_edge - (bags - 1) * batch
        bottom = np.array(margin < p1).astype(int)
        top = np.array(margin > batch - p1).astype(int)
        costs_new = costs * ((1 - np.exp((margin - p1) * bottom)) + 
        (1 - np.exp((batch - margin - p1) * top)))
    return costs_new

def iter_trans(costs, surp, lack, batch = 20, stair = 0):
    N = len(costs) # number of total area
    K = len(surp[0]) # type number of total cloth

    # the solution of last iteration
    X_last = np.zeros((N,N,K))
    costs_update = costs / batch
    obj_last = np.sum(costs_edge)

    iter_num = 0

    while 1:
        '''
        surp_inc = surp - np.sum(X_last, axis=1).reshape(N,K)
        lack_inc = lack - np.sum(X_last, axis=0).reshape(N,K)
        # solve each subproblem for every type
        res = muti_trans(costs_update, surp_inc, lack_inc)
        X_inc = res['var']
        X = X_inc + X_last
        '''
        res = muti_trans(costs_update, surp, lack)
        X = res['var']
        
        # total traffic volume of all type of goods
        # in each edge
        x_edge = np.sum(X, axis=2)
        obj = np.sum(np.ceil(x_edge / batch) * costs)

        if abs(obj - obj_last) < 0.0001 or iter_num > 100:
            print(obj - obj_last)
            break

        obj_last = obj
        X_last = X
        iter_num = iter_num + 1

        # update the costs of linear model
        costs_update = updata_costs(costs, batch, x_edge, stair)

    return X_last

def intp_trans(costs, surp, lack, batch = 20):
    N = len(costs)
    K = len(surp[0])

    surp_add = np.zeros((1,K))

    for k in range(K):
        # target nodes
        l_k = lack[:,k]
        
        # source nodes
        s_k = surp[:,k]

        T_s = np.sum(s_k) - np.sum(l_k)

        if T_s < 0: # the surplus is less than lack
            surp_add[0,k] = -T_s

    surplus = np.vstack((surplus, surp_add))

    Rows = range(N + 1)
    Cols = range(N)
    Types = range(K)

    prob = LpProblem("Transportation Problem as 0-1", LpMinimize)

    # trans_id = [[[(i,j,k) for k in range(K)] for j in range(N)] for i in range(N)]
    X_ijk = LpVariable.dicts('trans', (Rows,Cols,Types), 
                            lowBound = 0)

    # path_id = [[(i,j) for j in range(N)] for i in range(N)]
    Y_ij = LpVariable.dicts('path', (Rows,Cols), 
                            lowBound = 0,
                            cat = LpInteger)

    prob += lpSum([lpDot([Y_ij[i][j] for j in Cols], costs[i]) for i in Rows])
    
    for k in Types:
        for i in Rows:
            prob += lpSum([X_ijk[i][j][k] for j in Cols]) <= surplus[i][k]

        for j in Cols:
            prob += lpSum([X_ijk[i][j][k] for i in Rows]) == lack[j][k]

    for i in Rows:
        for j in Cols:
            prob += lpSum([X_ijk[i][j][k] for k in Types]) <= batch * Y_ij[i][j]
    
    prob.solve()

    obj = value(prob.objective)
    if LpStatus[prob.status] != 'Optimal':
        print("Status:", LpStatus[prob.status])
    X = np.array([[[value(X_ijk[i][j][k]) for k in Types] for j in Cols] for i in Rows])
    X = X[0:N,:,:]
    return X

def des_trans(costs, surp, lack, batch = 20):
    N = len(costs)
    K = len(surp[0])

    return np.zeros((N,N,K))