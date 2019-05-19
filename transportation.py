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

def bag_trans(surp, lack, bags):
    N = len(surp) # number of total area
    K = len(surp[0]) # type number of total cloth

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
    res = muti_trans(costs_update, surp, lack)
    X_last = res['var']
    '''
    obj_last = np.sum(costs_edge)

    iter_num = 0

    while 1:
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
    '''
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

def PSO_trans(costs, surp, lack, batch = 20):
    psize = 50 
    iters = 100
    w = 0.6
    vmax = 5
    c1 = 2
    c2 = 2

    N = len(costs) # number of total area
    K = len(surp[0]) # type number of total cloth

    costs_update = costs / batch
    res = muti_trans(costs_update, surp, lack)
    X_init = res['var']

    # initialize particles
    particles = []
    vels = []
    fits = np.zeros((psize,))
    n_id = [i for i in range(N)]
    for i in range(psize):
        vel = (np.random.rand(N,N,K) - 0.5) * 2 * v_max * (i > 0)
        '''
        x_rand = X_init
        penalty = 0
        for k in range(K):
            vel[n_id,n_id,k] = 0
            vel_s = np.sum(vel[:,:,k], axis=0)
            vel[:,:,k] -= np.repeat(vel_s.reshape(1,N), N, axis=0) / (N-1)
            vel[n_id,n_id,k] = 0
            x_rand[:,:,k] += vel[:,:,k]
            s_rand = surp[:,k] - np.sum(x_rand[:,:,k], axis=1)
            if np.min(s_rand) < 0 or np.min(x_rand[:,:,k]) < 0:
                penalty += 1000000
        '''
        vel[n_id,n_id,:] = 0
        vel_s = np.sum(vel[:,:,k],axis=0).reshape((1,N,K))
        vel -= np.repeat(vel_s, N, axis=0) / (N - 1)
        vel[n_id,n_id,:] = 0
        x_rand = X_init + vel
        s_rand = surp - np.sum(x_rand, axis=1)
        penalty = np.max(np.max(np.array(s_rand < 0).astype(int), axis=0)
         + np.max(np.array(x_rand < 0).astype(int), axis=(0,1))) * 100000.0
        particles.append(x_rand)
        vels.append(vel)
        y_rand = np.sum(x_rand, axis=2)
        fits[i] = np.sum(np.ceil(y_rand / batch) * costs) + penalty

    pbest_f = fits
    pbest = particles
    gbest_f = np.min(pbest_f)
    i = np.where(pbest_f == gbest_f)[0][0]
    gbest = particles[i]

    iter_num = 0
    while iter_num < iters:
        for i in range(psize):
            vels[i] = (w * vels[i] + c1 * random.random() * (pbest[i] - particles[i])
             + c2 * random.random() * (gbest - particles[i]))
            particles[i] += vels[i]
            xi = particles[i]

            s_rand = surp - np.sum(xi, axis=1)
            penalty = np.max(np.max(np.array(s_rand < 0).astype(int), axis=0)
            + np.max(np.array(xi < 0).astype(int), axis=(0,1))) * 100000.0
            '''
            penalty = 0
            for k in range(K):
                s_rand = surp[:,k] - np.sum(xi[:,:,k], axis=1)
                if np.min(s_rand) < 0 or np.min(xi[:,:,k]) < 0:
                    penalty += 1000000
            '''
            y_rand = np.sum(xi, axis=2)
            fits[i] = np.sum(np.ceil(y_rand / batch) * costs) + penalty
            if fits[i] < pbest_f[i]:
                pbest_f[i] = fits[i]
                pbest[i] = particles[i]
                if fits[i] < gbest_f:
                    gbest_f = fits[i]
                    gbest = particles[i]

    return gbest

def intra_alloc(sale_rate, storage, period = 2, min_depth = 3):
    # sale_rate: (N,); storage: (N,)
    N = len(sale_rate)
    sale_rate = sale_rate.reshape((N,))
    storage = sale_rate.reshape((N,))

    depth = storage / (sale_rate + 0.01 * np.array(sale_rate == 0).astype(int))
    shorts = min_depth * sale_rate - storage
    shorts = shorts * np.array(shorts > 0).astype(int)
    runds = 2 * period * sale_rate - storage
    runds = runds * np.array(runds > 0).astype(int)

    x_trans = muti_trans(np.zeros((N,N)), runds.reshape((N,1)), shorts.reshape((N,1)))
    x_trans = np.sum(x_trans, axis=2)
    # the (N+1)th means the warehouse
    result_x = np.zeros((N+1,N))
    result_x[0:N,:] = x_trans

    result_x[N,:] = (runds - np.sum(x_trans, axis=0)).reshape((1,N))

    return result_x


    
