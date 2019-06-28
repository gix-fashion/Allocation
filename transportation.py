# -*- coding: utf-8 -*-
from pulp import *
import numpy as np
import math
import random
import Queue

def pairwise(iterable):
    '''
    Given an iterable, yield the items in it in pairs. For instance:

        list(pairwise([1,2,3,4])) == [(1,2), (3,4)]
    '''
    x = iter(iterable)
    return zip(x, x)

class Simplex():

    def __init__(self, A = [], b = [], c = [], B = []):
        self._A = A # 系数矩阵
        if np.min(A) < 0:
            print('err in init')
            print(A)
            print('err' + i)
        self._b = b #
        self._c = c #
        self._B = B #基变量的下标集合
        self.row = len(b) #约束个数

    def solve(self):
        #假设线性规划形式是标准形式（都是等式）
        A = self._A
        b = self._b
        c = self._c
        B = self._B
        (x, B, obj, A_new, b_new) = self.Simplex(A, b, c)

        err = np.dot(self._A, x) - self._b
        ccc = True
        if np.sum(err * err) > 1e-6:
            ccc = False

        return (x, B, obj, ccc)

    def setProblem(self, A, b, c):
        self._A = A
        self._b = b
        self._c = c
        self.row = len(b)
    
    def setBase(self, B):
        self._B = B

    #算法入口
    def Simplex(self,A,b,c):
        B = self._B
        m_B = A[:,B]
        A_op = np.zeros((self.row, len(c)))
        for i in range(len(A[0])):
            A_op[:,i] = np.linalg.solve(m_B, A[:,i])
        b_op = np.linalg.solve(m_B, b) 

        #函数目标值
        obj = np.dot(c[B],b_op)

        c = np.dot(c[B].reshape((1,-1)),A_op) - c
        c = c[0]

        # entering basis
        e = np.argmax(c)
        # 找到最大的检验数，如果大于0，则目标函数可以优化
        while c[e] > 0:
            theta = []
            for i in range(len(b_op)):
                if A_op[i][e] > 0:
                    theta.append(b_op[i] / A_op[i][e])
                else:
                    theta.append(float("inf"))

            l = np.argmin(np.array(theta))

            if theta[l] == float('inf'):
                print 'unbounded'
                return False

            (B, A_op, b_op, c, obj) = self._PIVOT(B, A_op, b_op, c, obj, l, e)

            e = np.argmax(c)

        x = self._CalculateX(B,A_op,b_op,c)
        # err = np.dot(A[:,B], x[B]) - b

        self._B = B
        return (x, B, obj, A_op, b_op)

    #得到完整解
    def _CalculateX(self,B,A,b,c):

        x = np.zeros(len(self._c),dtype=float)
        x[B] = b
        return x

    # 基变换
    def _PIVOT(self,B,A,b,c,z,l,e):
        # main element is a_le
        # l represents leaving basis
        # e represents entering basis

        main_elem = A[l][e]
        #scaling the l-th line
        A[l] = A[l]/main_elem
        b[l] = b[l]/main_elem

        #change e-th column to unit array
        for i in range(self.row):
            if i != l:
                b[i] = b[i] - A[i][e] * b[l]
                A[i] = A[i] - A[i][e] * A[l]

        #update objective value
        z -= b[l]*c[e]

        c = c - c[e] * A[l]

        # change the basis
        B[l] = e

        return (B, A, b, c, z)

class TableWork():
    def __init__(self, A = [], b = [], c = [], B = []):
        self._A = A # 成本矩阵
        self._b = b # 产量
        self._c = c # 销量
        self._B = B # 基变量
        self._u = []
        self._v = []
        self.row = len(A)
        self.column = len(A[0])

    def initialize(self):
        As = np.zeros((self.row, self.column))
        bs = np.zeros((self.row,))
        bs[:] = self._b
        cs = np.zeros((self.column,))
        cs[:] = self._c
        self._B = np.zeros((self.row, self.column))
        x = np.zeros((self.row, self.column))
        while np.min(As) == 0:
            min_c = np.argmin(self._A + As).astype(int)
            new_B = [min_c / self.column, min_c % self.column]
            self._B[new_B[0]][new_B[1]] = 1
            if bs[new_B[0]] < cs[new_B[1]]:
                As[new_B[0],:] = 1e6
                cs[new_B[1]] -= bs[new_B[0]]
                x[new_B[0]][new_B[1]] = bs[new_B[0]]
                bs[new_B[0]] = 0
                # print(str(new_B[0]) + ', ' + str(new_B[1]) + ', b')
            elif bs[new_B[0]] > cs[new_B[1]]:
                As[:,new_B[1]] = 1e6
                bs[new_B[0]] -= cs[new_B[1]]
                x[new_B[0]][new_B[1]] = cs[new_B[1]]
                cs[new_B[1]] = 0
                # print(str(new_B[0]) + ', ' + str(new_B[1]) + ', c')
            
            else:
                x[new_B[0]][new_B[1]] = cs[new_B[1]]
                cs[new_B[1]] = 0
                bs[new_B[0]] = 0
                sup_pos = np.where((1 - self._B[new_B[0],:]) * (1 - As[new_B[0],:]) == 1)[0]
                if len(sup_pos) > 0:
                    self._B[new_B[0]][sup_pos[0]] = 1
                else:
                    sup_pos = np.where((1 - self._B[:,new_B[1]]) * (1 - As[:,new_B[1]]) == 1)[0]
                    if len(sup_pos) > 0:
                        self._B[sup_pos[0]][new_B[1]] = 1
                As[:,new_B[1]] = 1e6
                As[new_B[0],:] = 1e6
            
        '''
        new_B = [np.argmax(bs), np.argmax(cs)]
        self._B[new_B[0]][new_B[1]] = 1
        x[new_B[0]][new_B[1]] = cs[new_B[1]]
        '''
        return (x, bs, cs)
    
    def setBase(self, B):
        self._B = B

    def base2x(self):
        bs = np.zeros((self.row,))
        bs[:] = self._b
        cs = np.zeros((self.column,))
        cs[:] = self._c
        x = np.zeros((self.row, self.column))
        base_s = np.where(self._B == 1)
        for i in range(len(base_s[0])):
            pos = (base_s[0][i], base_s[1][i])
            if bs[pos[0]] <= cs[pos[1]]:
                x[pos[0]][pos[1]] = bs[pos[0]]
                cs[pos[1]] -= bs[pos[0]]
                bs[pos[0]] = 0
            else:
                bs[pos[0]] -= cs[pos[1]]
                x[pos[0]][pos[1]] = cs[pos[1]]
                cs[pos[1]] = 0
                
        return x

    def testNumber(self):
        W = np.zeros((self.row, self.column))
        W[:,:] = self._A
        B_c = np.zeros((self.row, self.column))

        u = np.zeros((self.row,))
        v = np.zeros((self.column))
        q = Queue.Queue(self.row + self.column)
        Base_in = np.where((self._B[0,:] == 1))[0]
        for ib in Base_in:
            if B_c[0][ib] == 0:
                q.put((0, ib, 0))
        
        k = 0
        while k < self.row + self.column - 1:
            head = q.get()
            if B_c[head[0]][head[1]] == 1:
                continue
            if head[2] == 0:
                v[head[1]] = self._A[head[0]][head[1]] - u[head[0]]
                W[:,head[1]] = W[:,head[1]] - v[head[1]]
                B_c[head[0]][head[1]] = 1
                in_Base = np.where((self._B[:,head[1]] == 1))[0]
                for ib in in_Base:
                    if B_c[ib][head[1]] == 0:
                        q.put((ib, head[1], 1))
                k = k + 1
            else:
                u[head[0]] = self._A[head[0]][head[1]] - v[head[1]]
                W[head[0],:] = W[head[0],:] - u[head[0]]
                B_c[head[0]][head[1]] = 1
                in_Base = np.where((self._B[head[0],:] == 1))[0]
                for ib in in_Base:
                    if B_c[head[0]][ib] == 0:
                        q.put((head[0], ib, 0))
                k = k + 1
        
        self._u = u
        self._v = v
        
        if np.max(np.fabs(np.repeat(u.reshape((self.row,1)), self.column, axis=1) + np.repeat(v.reshape((1,self.column)),
        self.row, axis=0) + W - self._A)) > 1e-8:
            print('testNum error')
            print(np.max(np.fabs(np.repeat(u.reshape((self.row,1)), self.column, axis=1) + np.repeat(v.reshape((1,self.column)),
            self.row, axis=0) + W - self._A)))
            # print('1' + W)
        
        if np.sum(W * self._B) != 0:
            print('testNum 0f base error')
            # print('1' + W)
        
        return W

    def adjust(self, pos, x):
        stack = []
        pos_ocp = np.ones((self.row, self.column))

        dd = 0
        cur_pos = pos
        
        # find cycle
        while 1:
            if dd == 0:
                if self._B[cur_pos[0]][pos[1]] == 1:
                    stack.append((cur_pos[0], pos[1], 0))
                    break
                nexts = np.where(self._B[cur_pos[0],:] == 1)[0]
                pp_n = []
                for n in nexts:
                    if pos_ocp[cur_pos[0]][n] == 0:
                        continue
                    if np.sum(self._B[:,n]) > 1:
                        pp_n.append((cur_pos[0], n, dd))
                        if x[cur_pos[0]][n] > 0:
                            stack.append((cur_pos[0], n, dd))
                            cur_pos = (cur_pos[0], n, dd)
                            pos_ocp[cur_pos[0]][n] = 0
                            dd = 1
                            break
                if dd == 0:
                    if len(pp_n) > 0:
                        rn = random.randint(0, len(pp_n) - 1)
                        stack.append(pp_n[rn])
                        cur_pos = pp_n[rn]
                        pos_ocp[cur_pos[0]][cur_pos[1]] = 0
                        dd = 1
                
                if dd == 0:
                    if len(stack) == 0:
                        print('cannot find a cycle  ' + str(dd))
                        print(pos)
                        print(self._B)
                        print(self._A)
                        print(x)
                        print('1' + x)
                        return False
                    stack.pop()
                    cur_pos = stack[-1]
                    dd = 1
            else:
                nexts = np.where(self._B[:, cur_pos[1]] == 1)[0]
                for n in nexts:
                    if pos_ocp[n][cur_pos[1]] == 0:
                        continue
                    if np.sum(self._B[n,:]) > 1:
                        stack.append((n, cur_pos[1], dd))
                        cur_pos = (n, cur_pos[1], dd)
                        pos_ocp[n][cur_pos[1]] = 0
                        dd = 0
                        break
                
                if dd == 1:
                    stack.pop()
                    if len(stack) == 0:
                        cur_pos = pos
                    else:
                        cur_pos = stack[-1]
                    dd = 0
            
            # print(cur_pos)

        theta = 1e8
        # str_out = '(' + str(pos[0]) + ', ' + str(pos[1]) + ', ' + '0) '
        for ps in stack:
            # str_out += '(' + str(ps[0]) + ', ' + str(ps[1]) + ', ' + str(ps[2]) + ',' + str(
            # x[ps[0]][ps[1]]) + ') '
            if ps[2] == 0:
                if x[ps[0]][ps[1]] < theta:
                    theta = x[ps[0]][ps[1]]
                    out_base = (ps[0], ps[1])
                elif x[ps[0]][ps[1]] == theta:
                    if random.random() > 0.5:
                        out_base = (ps[0], ps[1])

        x[pos[0]][pos[1]] = theta
        # str2 = str(theta) + 'after: (' + str(pos[0]) + ', ' + str(pos[1]) + ', ' + str(
        #     x[pos[0]][pos[1]]
        # )
        for ps in stack:
            x[ps[0]][ps[1]] += theta * (ps[2] * 2 - 1)
            # str2 += '(' + str(ps[0]) + ', ' + str(ps[1]) + ', ' + str(
            #     x[ps[0]][ps[1]]
            # )
        
        self._B[pos[0]][pos[1]] = 1
        self._B[out_base[0]][out_base[1]] = 0

        # print(str_out)
        # print(str2)

        return x

    def solve(self, init_s = 1, x_init = []):
        if init_s == 1:
            (x, bs, cs) = self.initialize()
        else:
            x = x_init
            
        if np.sum(self._B) != self.row + self.column - 1 or np.min(np.sum(self._B,axis=0)
        ) < 1 or np.min(np.sum(self._B,axis=1)) < 1:
            print('base num: ' + str(np.sum(self._B)) + ', size:' + str(self.row)
            + ', ' + str(self.column))
            print(self._A)
            print(self._B)
            print(self._b)
            print(self._c)
            print(bs)
            print(cs)
            print(x)
            print('1' + x)
        
        while 1:
            W = self.testNumber()
            # print('test')
            min_w = np.argmin(W).astype(int)
            pos = [min_w / self.column, min_w % self.column]
            # print(W[pos[0]][pos[1]])
            if W[pos[0]][pos[1]] > -1e-10:
                break
            
            x = self.adjust(pos, x)
            # print('adjust')

        sum_raw = np.sum(x, axis=1)
        sum_column = np.sum(x, axis=0)
        if np.sum((sum_raw - self._b) ** 2) > 1e-8 or np.sum((sum_column - self._c) ** 2) > 1e-8:
            # print(self._b)
            # print(self._c)
            # print(x)
            print('constrait ee')
            print('1' + x)

        if np.min(x) < 0:
            print('min x neg')
            print('1' + x)
        return (x, self._B)

def trans(costs, x_max, y_max):
    row = len(costs)
    if row == 0:
        return {'objective': 0, 'var': []}
    col = len(costs[0])
    Rows = range(row)
    Cols = range(col)

    prob = LpProblem("Transportation Problem", LpMinimize)
    vars = LpVariable.dicts("x", (Rows, Cols), lowBound=0)

    # prob += lpSum([lpSum([vars[i][j] * costs[i][j] for j in Cols]) for i in Rows])
    prob += lpSum([lpDot([vars[i][j] for j in Cols], costs[i]) for i in Rows])

    for i in range(row):
        prob += lpSum([vars[i][j] for j in Cols]) == x_max[i]

    for j in range(col):
        prob += lpSum([vars[i][j] for i in Rows]) == y_max[j]

    prob.solve()
    if LpStatus[prob.status] != 'Optimal':
        print("Status:", LpStatus[prob.status])
    return {'objective': value(prob.objective), 'var': [[value(vars[i][j]) for j in range(col)] for i in range(row)]}

def muti_trans(costs, surp, lack):
    N = len(costs)  # number of total area
    K = len(surp[0])  # type number of total cloth

    X = np.zeros((N, N, K))
    obj = 0
    for i in range(K):
        # target nodes
        l_k = lack[:, i]

        # source nodes
        s_k = surp[:, i]

        T_s = np.sum(s_k) - np.sum(l_k)

        index_l = np.where(l_k > 0)[0]
        index_s = np.where(s_k > 0)[0]
        if len(index_l) == 0 or len(index_s) == 0:
            continue

        # costs of k type
        costs_k = costs[index_s, :]
        costs_k = costs_k[:, index_l]
        costs_k = np.reshape(costs_k, (len(index_s), len(index_l)))
        # lack of k type
        lack_k = l_k[index_l]
        surplus_K = s_k[index_s]

        if T_s == 0:
            # solve transportation problem
            res = trans(costs_k, surplus_K, lack_k)
            x_k = np.array(res['var'])
            obj_k = res['objective']
        elif T_s < 0:  # the surplus is less than lack
            cost_add = np.zeros((1, len(index_l)))
            costs_k = np.vstack((costs_k, cost_add))
            surplus_K = np.hstack((surplus_K, -T_s))
            res = trans(costs_k, surplus_K, lack_k)
            x_k = np.array(res['var'])
            x_k = x_k[0:len(index_s), :]
            obj_k = res['objective']
        else:
            cost_add = np.zeros((len(index_s), 1))
            costs_k = np.hstack((costs_k, cost_add))
            lack_k = np.hstack((lack_k, T_s))
            res = trans(costs_k, surplus_K, lack_k)
            x_k = np.array(res['var'])
            x_k = x_k[:, 0:len(index_l)]
            obj_k = res['objective']

        for li in range(len(index_l)):
            X[index_s, index_l[li], i] = x_k[:, li]
        if obj_k != None:
            obj += obj_k

    return {'objective': obj, 'var': X}

def muti_trans2(costs, costs_SKC, lack_SKC,
    surp_SKC, type_SKC, LI_SKC, SI_SKC, N, K):
    
    X = np.zeros((N, N, K))
    obj = 0
    x_SKC = []
    for i in range(len(costs_SKC)):
        # costs of k type
        costs_k = costs_SKC[i]
        index_s = SI_SKC[i]
        index_l = LI_SKC[i]
        costs_n = costs[index_s, :]
        costs_n = costs_n[:, index_l]
        costs_k[0:len(index_s),0:len(index_l)] = costs_n

        res = trans(costs_k, surp_SKC[i], lack_SKC[i])
        x_k = np.array(res['var'])
        x_SKC.append(x_k)
        x_k = x_k[0:len(index_s), 0:len(index_l)]
        obj_k = res['objective']

        for li in range(len(index_l)):
            X[index_s, index_l[li], type_SKC[i]] = x_k[:, li]
        if obj_k != None:
            obj += obj_k

    return {'objective': obj, 'var': X, 'var_SKC': x_SKC}

def pretreat(costs_SKC, lack_SKC, surp_SKC, x_SKC):
    As = []
    bs = []
    Bs = []
    for i in range(len(costs_SKC)):
        # costs of k type
        costs_k = costs_SKC[i]

        ns = len(costs_k)
        nl = len(costs_k[0])

        lp_A1 = np.repeat(np.eye(ns), nl, axis=1)
        lp_A2 = np.tile(np.eye(nl), (1,ns))
        lp_A = np.vstack((lp_A1, lp_A2))
        lp_A = lp_A[0:(ns+nl-1),:] 
        As.append(lp_A)

        lp_b = np.hstack((surp_SKC[i],lack_SKC[i]))
        bs.append(lp_b[0:(ns+nl-1)])

        xi = x_SKC[i].reshape((1,-1))
        lp_B = np.where(xi > 0)[1]
        if len(lp_B) < ns + nl - 1:
            rn = ns+nl-len(lp_B)-1
            k = 0
            lp_N = np.where(xi == 0)[1]
            for mm in range(rn):
                while 1:
                    lp_B_test = np.hstack((lp_B, lp_N[k]))
                    if np.linalg.matrix_rank(lp_A[:,lp_B_test]) == len(lp_B) + 1:
                        lp_B = lp_B_test
                        k += 1
                        break
                    k += 1
        Bs.append(lp_B)

    return (As, bs, Bs)

def pretreat2(costs, surp, lack):
    costs_SKC = []
    lack_SKC = []
    surp_SKC = []
    type_SKC = []
    LI_SKC = []
    SI_SKC = []
    K = len(surp[0])
    for i in range(K):
        # target nodes
        l_k = lack[:, i]

        # source nodes
        s_k = surp[:, i]

        T_s = np.sum(s_k) - np.sum(l_k)

        index_l = np.where(l_k > 0)[0]
        index_s = np.where(s_k > 0)[0]
        if len(index_l) == 0 or len(index_s) == 0:
            continue

        type_SKC.append(i)
        LI_SKC.append(index_l)
        SI_SKC.append(index_s)
        # costs of k type
        costs_k = costs[index_s, :]
        costs_k = costs_k[:, index_l]
        costs_k = np.reshape(costs_k, (len(index_s), len(index_l)))
        # lack of k type
        lack_k = l_k[index_l]
        surplus_K = s_k[index_s]

        if T_s < 0:  # the surplus is less than lack
            cost_add = np.zeros((1, len(index_l)))
            costs_k = np.vstack((costs_k, cost_add))
            surplus_K = np.hstack((surplus_K, -T_s))
        elif T_s > 0:
            cost_add = np.zeros((len(index_s), 1))
            costs_k = np.hstack((costs_k, cost_add))
            lack_k = np.hstack((lack_k, T_s))
        
        costs_SKC.append(costs_k)
        lack_SKC.append(lack_k)
        surp_SKC.append(surplus_K)
    return (costs_SKC, lack_SKC, surp_SKC, type_SKC, LI_SKC, SI_SKC)

def Simplex_trans(costs, costs_SKC, As, bs,
    type_SKC, LI_SKC, SI_SKC, Bs, N, K):
    
    X = np.zeros((N, N, K))
    obj = 0
    for i in range(len(costs_SKC)):
        # costs of k type
        costs_k = costs_SKC[i]
        index_s = SI_SKC[i]
        index_l = LI_SKC[i]
        costs_n = costs[index_s, :]
        costs_n = costs_n[:, index_l]
        costs_k[0:len(index_s),0:len(index_l)] = costs_n

        ns = len(costs_k)
        nl = len(costs_k[0])

        lp_A = As[i]

        lp_b = bs[i]

        lp_c = costs_k.reshape((1,-1))
        lp_c = lp_c[0]

        lp_B = Bs[i]

        S = Simplex(lp_A, lp_b, lp_c, lp_B)
        (xi, lp_B, obj_k, ccc) = S.solve()
        # Bs[i] = lp_B
        
        x_k = xi.reshape((ns,nl))
        x_k = x_k[0:len(index_s), 0:len(index_l)]

        for li in range(len(index_l)):
            X[index_s, index_l[li], type_SKC[i]] = x_k[:, li]
        if obj_k != None:
            obj += obj_k

    return {'objective': obj, 'var': X, 'var_Bs': Bs}

def bag_trans(surp, lack, bags):
    N = len(surp)  # number of total area
    K = len(surp[0])  # type number of total cloth

    Rows = range(N)
    Cols = range(N)
    Types = range(K)

    prob = LpProblem("Transportation Problem as 0-1", LpMinimize)

    # trans_id = [[[(i,j,k) for k in range(K)] for j in range(N)] for i in range(N)]
    X_ijk = LpVariable.dicts('trans', (Rows, Cols, Types),
                             lowBound=0)

    prob += X_ijk[0][0][0]

    for k in Types:
        for i in Rows:
            prob += lpSum([X_ijk[i][j][k] for j in Cols]) <= surp[i][k]

        for j in Cols:
            prob += lpSum([X_ijk[i][j][k] for i in Rows]) == lack[j][k]

    for i in Rows:
        for j in Cols:
            prob += lpSum([X_ijk[i][j][k] for k in Types]) <= bags[i][j]

    prob.solve()

    if LpStatus[prob.status] != 'Optimal':
        # print("Status:", LpStatus[prob.status])
        return (False, 0)
    else:
        # print("valid")
        return (True, [[[X_ijk[i][j][k] for k in Types] for j in Cols] for i in Rows])

def updata_costs(costs, batch, x_edge, stair, p1=3, p2=1):
    if stair == 0:
        bags = np.ceil(x_edge / batch)
        margin = x_edge - (bags - 1) * batch
        bottom = np.array(margin < p1).astype(int)
        top = np.array(margin > batch - p1).astype(int)
        costs_new = costs * ((1 - np.exp((margin - p1) * bottom)) +
                             (1 - np.exp((batch - margin - p1) * top)))
    elif stair == 1:
        pen_edge = np.array(x_edge < p1).astype(int)
        pen_dom = x_edge * (1 - pen_edge) + (x_edge * 0.1 + 0.001) * pen_edge
        costs_new = costs / pen_dom
    return costs_new

def iter_trans(costs, surp, lack, batch=200, stair=1):
    N = len(costs)  # number of total area
    K = len(surp[0])  # type number of total cloth

    # the solution of last iteration
    X_last = np.zeros((N, N, K))
    costs_update = costs
    obj_last = np.sum(costs)

    (costs_SKC, lack_SKC, surp_SKC, type_SKC, LI_SKC, SI_SKC) = pretreat2(costs, surp, lack)
    
    iter_num = 0

    while 1:
        # res = muti_trans(costs_update, surp, lack)
        res = muti_trans2(costs_update, costs_SKC, lack_SKC, surp_SKC,
        type_SKC, LI_SKC, SI_SKC, N, K)
        X = res['var']

        # total traffic volume of all type of goods
        # in each edge
        x_edge = np.sum(X, axis=2)
        obj = np.sum(np.array(x_edge > 0).astype(int) * costs)
        bags = np.sum(np.array(x_edge > 0).astype(int))
        print('obj:' + str(obj) + ", bags: " + str(bags) + ", average: "
              + str(np.sum(x_edge) / bags))

        if abs(obj - obj_last) < 0.0001 or iter_num > 100:
            print(obj - obj_last)
            break

        obj_last = obj
        X_last = X
        iter_num = iter_num + 1

        # update the costs of linear model
        costs_update = updata_costs(costs, batch, x_edge, stair)

    obj_best = obj
    print('obj_best: ' + str(obj_best))
    
    iter_num = 0
    anni_state = np.zeros((N, N))
    x_SKC = res['var_SKC']
    (As, bs, Bs) = pretreat(costs_SKC, lack_SKC, surp_SKC, x_SKC)
    
    while iter_num < 100:
        costs_update = costs_update * np.array(x_edge > 0).astype(int) * (1 - anni_state)
        max_cost = np.max(costs_update)
        max_edge = np.where(costs_update == max_cost)
        anni_state[max_edge[0][0]][max_edge[1][0]] = 1

        costs_update = updata_costs(costs, batch, x_edge, stair)
        costs_update = costs_update * (1 - anni_state)

        res = Simplex_trans(costs_update, costs_SKC, As, bs, type_SKC, LI_SKC, SI_SKC, Bs, N, K)

        # res = muti_trans(costs_update, surp, lack)
        # res = muti_trans2(costs_update, costs_SKC, lack_SKC, surp_SKC,
        # type_SKC, LI_SKC, SI_SKC, N, K)
        X = res['var']
        x_edge = np.sum(X, axis=2)
        # Bs = res['var_Bs']

        obj = np.sum(np.array(x_edge > 0).astype(int) * costs)
        bags = np.sum(np.array(x_edge > 0).astype(int))
        print('obj:' + str(obj) + ", bags: " + str(bags) + ", average: "
              + str(np.sum(x_edge) / bags))

        if obj < obj_best:
            obj_best = obj
            print('obj_best: ' + str(obj_best))

        iter_num += 1
        if np.sum(np.array(x_edge > 0).astype(int) * (1 - anni_state)) == 0:
            break

    s2s_skc = np.sum(np.array(X > 0).astype(int), axis=1)
    invalid_num = np.sum(np.array(s2s_skc > 2).astype(int))
    print('invalid_num = ' + str(invalid_num) + ',  max_: ' + str(np.max(s2s_skc)))
    small_bags = np.sum(np.array(x_edge < 2).astype(int) * np.array(x_edge > 0).astype(int))
    print('small bag: ' + str(small_bags))
    
    return X

def Table_trans(costs, costs_SKC, lack_SKC,
    surp_SKC, type_SKC, LI_SKC, SI_SKC, Bs_last, x_last, N, K, setB=1):
    X = np.zeros((N, N, K))
    obj = 0
    if setB == 0:
        Bs_last = []
        x_last = []
    for i in range(len(costs_SKC)):
    # for i in [73]:
        # costs of k type
        costs_k = costs_SKC[i]
        index_s = SI_SKC[i]
        index_l = LI_SKC[i]
        costs_n = costs[index_s, :]
        costs_n = costs_n[:, index_l]
        costs_k[0:len(index_s),0:len(index_l)] = costs_n

        Tw = TableWork(costs_k, surp_SKC[i], lack_SKC[i])
        
        if setB == 1:
            Tw.setBase(Bs_last[i])
            (x_k, B_k) = Tw.solve(init_s=0,x_init=x_last[i])
            Bs_last[i] = B_k
            x_last[i] = x_k
            # Tw.setBase(Bs_last[0])
            # (x_k, B_k) = Tw.solve(init_s=0,x_init=x_last[0])
        else:
            (x_k, B_k) = Tw.solve()
            Bs_last.append(B_k)
            x_last.append(x_k)
        
        # print(str(i) + ', size: ' + str(len(index_s)) + '*' + str(len(index_l)))

        # print('1' + x_k)

        x_k = x_k[0:len(index_s), 0:len(index_l)]

        for li in range(len(index_l)):
            X[index_s, index_l[li], type_SKC[i]] = x_k[:, li]

    return (X, Bs_last, x_last)

def iter_table(costs, surp, lack, max_iter=100):
    N = len(costs)  # number of total area
    K = len(surp[0])  # type number of total cloth

    # the solution of last iteration
    X_last = np.zeros((N, N, K))
    costs_update = costs
    obj_last = np.sum(costs)

    (costs_SKC, lack_SKC, surp_SKC, type_SKC, LI_SKC, SI_SKC) = pretreat2(costs, surp, lack)
    
    iter_num = 0
    setB = 0
    Bs = []
    x_SKC = []

    while 1:
        # res = muti_trans(costs_update, surp, lack)
        (X, Bs, x_SKC) = Table_trans(costs_update, costs_SKC, lack_SKC, surp_SKC,
        type_SKC, LI_SKC, SI_SKC, Bs, x_SKC, N, K, setB)

        setB = 1

        # total traffic volume of all type of goods
        # in each edge
        x_edge = np.sum(X, axis=2)
        obj = np.sum(np.array(x_edge > 0).astype(int) * costs)
        bags = np.sum(np.array(x_edge > 0).astype(int))
        print('obj:' + str(obj) + ", bags: " + str(bags) + ", average: "
              + str(np.sum(x_edge) / bags))

        if abs(obj - obj_last) < 0.0001 or iter_num > 100:
            print(obj - obj_last)
            break

        obj_last = obj
        X_last = X
        iter_num = iter_num + 1

        # update the costs of linear model
        costs_update = updata_costs(costs, 200, x_edge, 1)

    obj_best = obj
    print('obj_best: ' + str(obj_best))
    
    iter_num = 0
    anni_state = np.zeros((N, N))
    
    while iter_num < max_iter:
        costs_update = costs_update * np.array(x_edge > 0).astype(int) * (1 - anni_state)
        max_cost = np.max(costs_update)
        max_edge = np.where(costs_update == max_cost)
        anni_state[max_edge[0][0]][max_edge[1][0]] = 1

        costs_update = updata_costs(costs, 200, x_edge, 1)
        costs_update = costs_update * (1 - anni_state)

        (X, Bs, x_SKC) = Table_trans(costs_update, costs_SKC, lack_SKC, surp_SKC,
        type_SKC, LI_SKC, SI_SKC, Bs, x_SKC, N, K, setB)
        
        x_edge = np.sum(X, axis=2)
        # Bs = res['var_Bs']

        obj = np.sum(np.array(x_edge > 0).astype(int) * costs)
        bags = np.sum(np.array(x_edge > 0).astype(int))
        print('obj:' + str(obj) + ", bags: " + str(bags) + ", average: "
              + str(np.sum(x_edge) / bags))

        if obj < obj_best:
            obj_best = obj
            print('obj_best: ' + str(obj_best))

        iter_num += 1
        if np.sum(np.array(x_edge > 0).astype(int) * (1 - anni_state)) == 0:
            break

    s2s_skc = np.sum(np.array(X > 0).astype(int), axis=1)
    invalid_num = np.sum(np.array(s2s_skc > 2).astype(int))
    print('invalid_num = ' + str(invalid_num) + ',  max_: ' + str(np.max(s2s_skc)))
    small_bags = np.sum(np.array(x_edge < 2).astype(int) * np.array(x_edge > 0).astype(int))
    print('small bag: ' + str(small_bags))
    
    return X

def IL_trans(costs, surp, lack, batch=20):
    N = len(costs)  # number of total area
    K = len(surp[0])  # type number of total cloth

    # the solution of last iteration
    costs_update = costs
    obj_last = np.sum(costs)

    iter_num = 0

    while 1:
        res = muti_trans(costs_update, surp, lack)
        X = res['var']

        # total traffic volume of all type of goods
        # in each edge
        x_edge = np.sum(X, axis=2)
        obj = np.sum(np.array(x_edge > 0).astype(int) * costs)
        bags = np.sum(np.array(x_edge > 0).astype(int))
        print('obj:' + str(obj) + ", bags: " + str(bags) + ", average: "
              + str(np.sum(x_edge) / bags))

        if abs(obj - obj_last) < 0.0001 or iter_num > 100:
            print(obj - obj_last)
            break

        obj_last = obj
        iter_num = iter_num + 1

        # update the costs of linear model
        costs_update = updata_costs(costs, batch, x_edge, stair)

    obj_best = obj
    print('obj_best: ' + str(obj_best))

    iter_num = 0
    alpha = 5
    beta = 0.9
    base_p = np.ones((N, N))
    anni_prob = np.zeros((N, N))
    anni_state = np.zeros((N, N))
    anni_base = 0.1
    gene_base = 0.1
    while iter_num < 100:
        rd_num = np.random.rand(N, N)
        anni_state += np.array(anni_state < 0.5).astype(int) * np.array(
            rd_num < anni_base * base_p
        ).astype(int)
        anni_state -= np.array(anni_state > 0.5).astype(int) * np.array(
            rd_num < gene_base / (base_p + anni_prob)
        ).astype(int)

        costs_update = updata_costs(costs, batch, x_edge, stair)
        costs_update = costs_update * (1 - anni_state) + anni_state * 1e8

        res = muti_trans(costs_update, surp, lack)
        X = res['var']
        x_edge = np.sum(X, axis=2)

        obj = np.sum(np.array(x_edge > 0).astype(int) * costs)
        bags = np.sum(np.array(x_edge > 0).astype(int))
        print('obj:' + str(obj) + ", bags: " + str(bags) + ", average: "
              + str(np.sum(x_edge) / bags))

        if obj < obj_last:
            anni_prob += anni_state * (obj_last - obj)
            anni_prob -= (1 - anni_state) * anni_prob * (1 - beta)
            if obj < obj_best:
                base_p += anni_state
                obj_best = obj
                print('obj_best: ' + str(obj_best))
        else:
            anni_prob = anni_prob * beta

        obj_last = obj
        iter_num += 1

    s2s_skc = np.sum(np.array(X > 0).astype(int), axis=1)
    invalid_num = np.sum(np.array(s2s_skc > 2).astype(int))
    print('invalid_num = ' + str(invalid_num) + ',  max_: ' + str(np.max(s2s_skc)))
    small_bags = np.sum(np.array(x_edge < 3).astype(int) * np.array(x_edge > 0).astype(int))
    print('small bag: ' + str(small_bags))
    return X


def intp_trans(costs, surp, lack, batch=20):
    N = len(costs)
    K = len(surp[0])

    surp_add = np.zeros((1, K))

    for k in range(K):
        # target nodes
        l_k = lack[:, k]

        # source nodes
        s_k = surp[:, k]

        T_s = np.sum(s_k) - np.sum(l_k)

        if T_s < 0:  # the surplus is less than lack
            surp_add[0, k] = -T_s

    surplus = np.vstack((surp, surp_add))

    Rows = range(N + 1)
    Cols = range(N)
    Types = range(K)

    prob = LpProblem("Transportation Problem as 0-1", LpMinimize)

    # trans_id = [[[(i,j,k) for k in range(K)] for j in range(N)] for i in range(N)]
    X_ijk = LpVariable.dicts('trans', (Rows, Cols, Types),
                             lowBound=0)

    # path_id = [[(i,j) for j in range(N)] for i in range(N)]
    Y_ij = LpVariable.dicts('path', (Rows, Cols),
                            lowBound=0,
                            cat=LpInteger)

    prob += lpSum([lpDot([Y_ij[i][j] for j in Cols], costs[i]) for i in Cols])

    for k in Types:
        for i in Rows:
            prob += lpSum([X_ijk[i][j][k] for j in Cols]) <= surplus[i][k]

        for j in Cols:
            prob += lpSum([X_ijk[i][j][k] for i in Rows]) == lack[j][k]

    for i in Rows:
        for j in Cols:
            prob += lpSum([X_ijk[i][j][k] for k in Types]) <= batch * Y_ij[i][j]

    prob.solve(GLPK())

    obj = value(prob.objective)
    print(obj)
    if LpStatus[prob.status] != 'Optimal':
        print("Status:", LpStatus[prob.status])
    X = np.array([[[value(X_ijk[i][j][k]) for k in Types] for j in Cols] for i in Rows])
    X = X[0:N, :, :]
    return X


def des_trans(costs, surp, lack, batch=20):
    N = len(costs)
    K = len(surp[0])

    return np.zeros((N, N, K))


def PSO_trans(costs, surp, lack, batch=20):
    psize = 5000  # particle swarm size
    iters = 500
    w = 1
    vmax = 0.5
    c1 = 2
    c2 = 2
    pc1 = 0.1
    pc2 = 0.1

    N = len(costs)  # number of total area
    K = len(surp[0])  # type number of total cloth

    costs_update = costs / batch
    res = muti_trans(costs_update, surp, lack)
    X_init = res['var']
    y_init = np.sum(X_init, axis=2)
    s_init = surp - np.sum(X_init, axis=1)
    l_init = lack - np.sum(X_init, axis=0)
    penalty = np.sum(pc1 * np.array(s_init < 0).astype(int) * np.exp(-s_init) +
                     pc2 * np.exp(np.abs(l_init)))
    fit = np.sum(np.array(y_init > 0).astype(int) * costs) + penalty

    # initialize particles
    particles = []
    particles.append((fit, X_init, np.zeros((N, N, K))))
    n_id = [i for i in range(N)]

    for i in range(psize - 1):
        vel = (np.random.rand(N, N, K) - 0.5) * 2 * vmax
        x_rand = X_init + (np.random.rand(N, N, K) - 0.5) * 2 * 2
        x_rand[n_id, n_id, :] = 0
        x_rand = x_rand * np.array(x_rand > 0).astype(int)
        s_rand = surp - np.sum(x_rand, axis=1)
        l_rand = lack - np.sum(x_rand, axis=0)
        penalty = np.sum(pc1 * np.array(s_rand < 0).astype(int) * np.exp(-s_rand) +
                         pc2 * np.exp(np.abs(l_rand)))

        y_rand = np.sum(x_rand, axis=2)
        fit = np.sum(np.array(y_rand > 0).astype(int) * costs) + penalty
        particles.append((fit, x_rand, vel))

    pbest = particles
    gbest = min(particles)
    print(gbest[0])

    iter_num = 0
    while iter_num < iters:
        iter_num += 1
        for i in range(psize):
            vel = particles[i][2]
            xi = particles[i][1]
            vel = (w * vel + c1 * random.random() * (pbest[i][1] - xi)
                   + c2 * random.random() * (gbest[1] - xi))
            xi += vel
            xi[n_id, n_id, :] = 0
            xi = xi * np.array(xi > 0).astype(int)

            s_rand = surp - np.sum(xi, axis=1)
            l_rand = lack - np.sum(xi, axis=0)
            penalty = np.sum(pc1 * np.array(s_rand < 0).astype(int) * np.exp(-s_rand) +
                             pc2 * np.exp(np.abs(l_rand)))

            y_rand = np.sum(xi, axis=2)
            fit = np.sum(np.array(y_rand > 0).astype(int) * costs) + penalty
            particles[i] = (fit, xi, vel)
            if fit < pbest[i][0]:
                print(str(pbest[i][0]) + '->' + str(fit))
                pbest[i] = particles[i]
                if fit < gbest[0]:
                    gbest = particles[i]
                    print('gbest: ' + str(gbest[0]))

    return gbest[1]


def GA_parents(dnas, parents_num, tournament_size=5):
    for i in range(parents_num):
        yield min(random.sample(dnas, tournament_size))


def GA_selection(dnas, parents_num):
    psize = len(dnas)
    for i in range(parents_num):
        yield dnas[min(int(-10 * math.log(1 - random.random())), psize - 1)]


def GA_mute(dna_array, N, mrate):
    mute_rnum = np.random.rand(N, N)
    result = dna_array - np.array(mute_rnum < 2.0 * mrate / 3.0).astype(int) + np.array(
        mute_rnum > 1.0 - mrate / 3.0
    ).astype(int)
    result = result * np.array(result > 0).astype(int)
    return result


def GA_trans(costs, surp, lack, batch=20, psize=50):
    # params of GA
    mutation_rate = 0.1
    elite_size = 2
    init_var = 2
    max_iters = 100

    N = len(costs)  # number of total area
    K = len(surp[0])  # type number of total cloth

    costs_update = costs / batch
    res = muti_trans(costs_update, surp, lack)
    X_init = res['var']

    # generate init DNAs
    dnas = []
    fits = []
    y_init = np.ceil(np.sum(X_init, axis=2) / batch)
    dnas.append(y_init)
    fits.append(np.sum(y_init * costs))
    print(np.sum(y_init * costs))
    total_edge = len(np.where(y_init > 0)[0])
    init_rate = total_edge / (N * N)

    n_id = [i for i in range(N)]
    array_pos = np.array([[i * N + j for j in range(N)] for i in range(N)])
    r_num = 0

    for i in range(psize - 1):
        delta = np.random.rand(N, N)
        dna = y_init + (np.array(delta < 0.1 * init_rate).astype(int)
                        - np.array(y_init > 0).astype(int) * np.array(delta < 0.05).astype(int))
        dna[n_id, n_id] = 0
        # test if the dna is valid
        bags = dna * batch
        is_valid = bag_trans(surp, lack, bags)[0]
        if is_valid:
            penalty = 0
            print(np.sum(dna * costs))
            r_num += 1
        else:
            penalty = np.sum(costs)

        dnas.append(dna)
        fits.append(np.sum(dna * costs) + penalty)

    best_fit = min(fits)
    print("best_fit: " + str(best_fit) + ", r_num: " + str(r_num))

    scored_dnas = [(fits[i], dnas[i]) for i in range(psize)]
    iter_num = 0
    while iter_num < max_iters:
        iter_num += 1
        r_num = 0
        sorted_dnas = sorted(scored_dnas, key=lambda sdna: sdna[0])
        # elite elements
        scored_dnas = sorted_dnas[:elite_size]
        print(len(sorted_dnas))
        for parent1, parent2 in pairwise(GA_selection(sorted_dnas, psize - elite_size)):
            # crossover parents
            point1, point2 = sorted(random.randint(0, N * N) for _ in range(2))
            mask = np.array(array_pos <= point1).astype(int) + np.array(
                array_pos > point2
            ).astype(int)
            child1 = GA_mute(parent1[1] * mask + parent2[1] * (1 - mask), N, mutation_rate)
            child2 = GA_mute(parent2[1] * mask + parent1[1] * (1 - mask), N, mutation_rate)
            child1[n_id, n_id] = 0
            is_valid = bag_trans(surp, lack, child1 * batch)[0]
            if is_valid:
                penalty = 0
                print(np.sum(child1 * costs))
                r_num += 1
            else:
                penalty = np.sum(costs)
            fit = np.sum(child1 * costs) + penalty
            scored_dnas.append((fit, child1))
            child2[n_id, n_id] = 0
            is_valid = bag_trans(surp, lack, child2 * batch)[0]
            if is_valid:
                penalty = 0
                print(np.sum(child2 * costs))
                r_num += 1
            else:
                penalty = 20 * np.sum(costs)
            fit = np.sum(child2 * costs) + penalty
            scored_dnas.append((fit, child2))

        best_dna = min(scored_dnas)
        print("best_fit: " + str(best_dna[0]) + ", second: " +
              str(scored_dnas[1][0]) + ", r_num: " + str(r_num))
    X = np.array(bag_trans(surp, lack, best_dna[1] * batch)[1])
    return X


def intra_alloc(sale_rate, storage, period=2, min_depth=3):
    # sale_rate: (N,); storage: (N,)
    N = len(sale_rate)
    sale_rate = sale_rate.reshape((N,))
    storage = sale_rate.reshape((N,))

    depth = storage / (sale_rate + 0.01 * np.array(sale_rate == 0).astype(int))
    shorts = min_depth * sale_rate - storage
    shorts = shorts * np.array(shorts > 0).astype(int)
    runds = 2 * period * sale_rate - storage
    runds = runds * np.array(runds > 0).astype(int)

    x_trans = muti_trans(np.zeros((N, N)), runds.reshape((N, 1)), shorts.reshape((N, 1)))
    x_trans = np.sum(x_trans, axis=2)
    # the (N+1)th means the warehouse
    result_x = np.zeros((N + 1, N))
    result_x[0:N, :] = x_trans

    result_x[N, :] = (runds - np.sum(x_trans, axis=0)).reshape((1, N))

    return result_x



