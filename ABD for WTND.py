from __future__ import division, print_function
import math
import gurobipy as GRBPY
import time
import pandas as pd

start = time.clock()


# awkward restriction for 'callback'
def cbwarehouse(model, where):
    if where == GRBPY.GRB.Callback.MIPSOL:
        if model._iter >= 1:
            for j in range(model._J):
                model._c23[j].rhs = sum(model._C1[l] * model.cbGetSolution(model._y1[j, l]) for l in range(model._L))

            for k in range(model._K):
                model._c24[k].rhs = sum(model._C2[l] * model.cbGetSolution(model._y2[k, l]) for l in range(model._L))

            for w in range(model._W):
                model._c25[w].rhs = sum(model._C3[l] * model.cbGetSolution(model._y3[w, l]) for l in range(model._L))

            for j in range(model._J):
                for l in range(model._L):
                    model._c411[j, l].rhs = - model._M * model.cbGetSolution(model._y1[j, l]) + model._M

            for j in range(model._J):
                for l in range(model._L):
                    model._c413[j, l].rhs = -model._M * model.cbGetSolution(model._y1[j, l])

            for k in range(model._K):
                for l in range(model._L):
                    model._c421[k, l].rhs = - model._M * model.cbGetSolution(model._y2[k, l]) + model._M

            for k in range(model._K):
                for l in range(model._L):
                    model._c423[k, l].rhs = -model._M * model.cbGetSolution(model._y2[k, l])

            for w in range(model._W):
                for l in range(model._L):
                    model._c431[w, l].rhs = - model._M * model.cbGetSolution(model._y3[w, l]) + model._M

            for w in range(model._W):
                for l in range(model._L):
                    model._c433[w, l].rhs = -model._M * model.cbGetSolution(model._y3[w, l])

            model._c5.rhs = sum(model._c1[j][l] * model.cbGetSolution(model._y1[j, l]) for j in range(model._J) for l in
                                range(model._L)) + \
                            sum(model._c2[k][l] * model.cbGetSolution(model._y2[k, l]) for k in range(model._K) for l in
                                range(model._L)) + \
                            sum(model._c3[w][l] * model.cbGetSolution(model._y3[w, l]) for w in range(model._W) for l in
                                range(model._L))

            for n in range(model._N):
                model._c51[n].rhs = sum(
                    model._cn1[n][j][l] * model.cbGetSolution(model._y1[j, l]) for j in range(model._J) for l in
                    range(model._L))

            for n in range(model._N):
                model._c52[n].rhs = sum(
                    model._cn2[n][k][l] * model.cbGetSolution(model._y2[k, l]) for k in range(model._K) for l in
                    range(model._L))

            for n in range(model._N):
                model._c53[n].rhs = sum(
                    model._cn3[n][w][l] * model.cbGetSolution(model._y3[w, l]) for w in range(model._W) for l in
                    range(model._L))

        model._sub.optimize()

        if model._sub.status == GRBPY.GRB.INFEASIBLE:
            print("Iteration: ", model._iter)
            print("Adding feasibility cut...\n")

            lazycut = sum(model._c21[i].farkasdual * model._D[i] for i in range(model._I)) + \
                      sum(model._c23[j].farkasdual * (sum(model._C1[l] * model._y1[j, l] for l in range(model._L))) for
                          j in range(model._J)) + \
                      sum(model._c24[k].farkasdual * (sum(model._C2[l] * model._y2[k, l] for l in range(model._L))) for
                          k in range(model._K)) + \
                      sum(model._c25[w].farkasdual * (sum(model._C3[l] * model._y3[w, l] for l in range(model._L))) for
                          w in range(model._W)) + \
                      sum(model._c34[i, j].farkasdual * model._tc1[i][j] for i in range(model._I) for j in
                          range(model._J)) + \
                      sum(model._c35[j, k].farkasdual * model._tc2[j][k] for j in range(model._J) for k in
                          range(model._K)) + \
                      sum(model._c36[j, w].farkasdual * model._tc3[j][w] for j in range(model._J) for w in
                          range(model._W)) + \
                      sum(model._c411[j, l].farkasdual * (- model._M * model._y1[j, l] + model._M) for j in
                          range(model._J) for l in range(model._L)) + \
                      sum(model._c413[j, l].farkasdual * (-model._M * model._y1[j, l]) for j in range(model._J) for l in
                          range(model._L)) + \
                      sum(model._c421[k, l].farkasdual * (- model._M * model._y2[k, l] + model._M) for k in
                          range(model._K) for l in range(model._L)) + \
                      sum(model._c423[k, l].farkasdual * (-model._M * model._y2[k, l]) for k in range(model._K) for l in
                          range(model._L)) + \
                      sum(model._c431[w, l].farkasdual * (- model._M * model._y3[w, l] + model._M) for w in
                          range(model._W) for l in range(model._L)) + \
                      sum(model._c433[w, l].farkasdual * (-model._M * model._y3[w, l]) for w in range(model._W) for l in
                          range(model._L)) + \
                      model._c5.farkasdual * (sum(
                model._c1[j][l] * model._y1[j, l] for j in range(model._J) for l in range(model._L)) + \
                                              sum(model._c2[k][l] * model._y2[k, l] for k in range(model._K) for l in
                                                  range(model._L)) + \
                                              sum(model._c3[w][l] * model._y3[w, l] for w in range(model._W) for l in
                                                  range(model._L))) + \
                      sum(model._c51[n].farkasdual * (
                          sum(model._cn1[n][j][l] * model._y1[j, l] for j in range(model._J) for l in range(model._L)))
                          for n in range(model._N)) + \
                      sum(model._c52[n].farkasdual * (sum(model._cn2[n][k][l] * model._y2[k, l] for k in range(model._K)
                                                          for l in range(model._L))) for n in range(model._N)) + \
                      sum(model._c53[n].farkasdual * (sum(model._cn3[n][w][l] * model._y3[w, l] for w in range(model._W)
                                                          for l in range(model._L))) for n in range(model._N)) + \
                      model._c6.farkasdual * model._C

            model.cbLazy(lazycut >= 0)

            model._iter += 1
        elif model._sub.status == GRBPY.GRB.OPTIMAL:
            if model._sub.objval > model.cbGetSolution(model._maxshipcost) + 1e-6:
                print("Iteration: ", model._iter)
                print("Adding optimality cut...\n")

                lazycut = sum(model._c21[i].pi * model._D[i] for i in range(model._I)) + \
                          sum(model._c23[j].pi * (sum(model._C1[l] * model._y1[j, l] for l in range(model._L))) for j in
                              range(model._J)) + \
                          sum(model._c24[k].pi * (sum(model._C2[l] * model._y2[k, l] for l in range(model._L))) for k in
                              range(model._K)) + \
                          sum(model._c25[w].pi * (sum(model._C3[l] * model._y3[w, l] for l in range(model._L))) for w in
                              range(model._W)) + \
                          sum(model._c34[i, j].pi * model._tc1[i][j] for i in range(model._I) for j in
                              range(model._J)) + \
                          sum(model._c35[j, k].pi * model._tc2[j][k] for j in range(model._J) for k in
                              range(model._K)) + \
                          sum(model._c36[j, w].pi * model._tc3[j][w] for j in range(model._J) for w in
                              range(model._W)) + \
                          sum(model._c411[j, l].pi * (- model._M * model._y1[j, l] + model._M) for j in range(model._J)
                              for l in range(model._L)) + \
                          sum(model._c413[j, l].pi * (-model._M * model._y1[j, l]) for j in range(model._J) for l in
                              range(model._L)) + \
                          sum(model._c421[k, l].pi * (- model._M * model._y2[k, l] + model._M) for k in range(model._K)
                              for l in range(model._L)) + \
                          sum(model._c423[k, l].pi * (-model._M * model._y2[k, l]) for k in range(model._K) for l in
                              range(model._L)) + \
                          sum(model._c431[w, l].pi * (- model._M * model._y3[w, l] + model._M) for w in range(model._W)
                              for l in range(model._L)) + \
                          sum(model._c433[w, l].pi * (-model._M * model._y3[w, l]) for w in range(model._W) for l in
                              range(model._L)) + \
                          model._c5.pi * (sum(
                    model._c1[j][l] * model._y1[j, l] for j in range(model._J) for l in range(model._L)) + \
                                          sum(model._c2[k][l] * model._y2[k, l] for k in range(model._K) for l in
                                              range(model._L)) + \
                                          sum(model._c3[w][l] * model._y3[w, l] for w in range(model._W) for l in
                                              range(model._L))) + \
                          sum(model._c51[n].pi * (sum(
                              model._cn1[n][j][l] * model._y1[j, l] for j in range(model._J) for l in range(model._L)))
                              for n in range(model._N)) + \
                          sum(model._c52[n].pi * (sum(model._cn2[n][k][l] * model._y2[k, l] for k in range(model._K)
                                                      for l in range(model._L))) for n in range(model._N)) + \
                          sum(model._c53[n].pi * (sum(model._cn3[n][w][l] * model._y3[w, l] for w in range(model._W)
                                                      for l in range(model._L))) for n in range(model._N)) + \
                          model._c6.pi * model._C

                model.cbLazy(model._maxshipcost >= lazycut)

                model._iter += 1
        else:
            model.terminate()


class WareHouse:
    def __init__(self):
        # initialize data
        self.I = 11
        self.J = 37
        self.K = 30
        self.W = 30
        self.L = 3
        self.N = 8
        self.M = 100000
        self.C1 = [50, 100, 300]
        self.C2 = [50, 100, 150]
        self.C3 = [50, 100, 150]
        self.Cbar1 = 24
        self.Cbar2 = 48
        self.d = 4.7755 * 1e-6
        self.d1 = [0.07, 0.14, 0.28]
        self.d2 = [5.95, 11.9, 17.85]
        self.d3 = [3.89, 7.78, 11.66]
        self.tc1 = []
        self.tc2 = []
        self.tc3 = []
        self.p1 = []
        self.p2 = []
        self.p3 = []
        self.P1 = []
        self.P2 = []
        self.P3 = []
        self.c1 = []
        self.c2 = []
        self.c3 = []
        self.cn1 = []
        self.cn2 = []
        self.cn3 = []
        self.s1 = [0.03, 0.05, 0.12]
        self.s2 = [0.01, 0.06, 0.1]
        self.s3 = [0.02, 0.05, 0.08]

        self.D = [458, 266, 485, 395, 192, 834, 230, 147, 169, 332, 501]  # 原始值

        self.sigma = ([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
        self.C = 2.5
        self.Gamma = 5.3

        # initialize variables and constraints
        self.y1 = {}
        self.y2 = {}
        self.y3 = {}
        self.f1 = {}
        self.f2 = {}
        self.f3 = {}
        self.x1 = {}
        self.x2 = {}
        self.x3 = {}
        self.alpha = {}
        self.beta = {}
        self.gamma = {}
        self.delt = {}
        self.epsilon = {}
        self.eta = {}
        self.theta = {}
        self.lamda = {}
        self.Omega1 = {}
        self.Omega2 = {}
        self.Omega3 = {}
        self.maxshipcost = {}

        self.c = {}
        self.u10 = {}
        self.v10 = {}
        self.u1 = {}
        self.v1 = {}
        self.u2 = {}
        self.v2 = {}
        self.u3 = {}
        self.v3 = {}
        self.beta1 = {}
        self.beta2 = {}
        self.beta3 = {}
        self.beta4 = {}
        self.beta5 = {}
        self.beta6 = {}

        self.c01 = {}
        self.c02 = {}

        self.c11 = {}
        self.c12 = {}
        self.c13 = {}
        self.c14 = {}

        self.c21 = {}
        self.c22 = {}
        self.c23 = {}
        self.c24 = {}
        self.c25 = {}
        self.c26 = {}
        self.c27 = {}
        self.c28 = {}

        self.c31 = {}
        self.c32 = {}
        self.c33 = {}
        self.c34 = {}
        self.c35 = {}
        self.c36 = {}

        self.c411 = {}
        self.c412 = {}
        self.c413 = {}
        self.c421 = {}
        self.c422 = {}
        self.c423 = {}
        self.c431 = {}
        self.c432 = {}
        self.c433 = {}
        self.c44 = {}

        self.c5 = {}
        self.c51 = {}
        self.c52 = {}
        self.c53 = {}
        self.c54 = {}
        self.c551 = {}
        self.c552 = {}
        self.c561 = {}
        self.c562 = {}
        self.c571 = {}
        self.c572 = {}
        self.c581 = {}
        self.c582 = {}
        self.c591 = {}
        self.c592 = {}
        self.c5101 = {}
        self.c5102 = {}

        self.c6 = {}

    def read(self, filename):
        # input data
        with open(filename, "r") as data:
            for i in range(self.I):
                column = data.readline().split()
                ltc1 = []
                for j in range(self.J):
                    ltc1.append(float(column[j]))
                self.tc1.append(ltc1)
                # print(d1)

            for j in range(self.J):
                column = data.readline().split()
                ltc2 = []
                for k in range(self.K):
                    ltc2.append(float(column[k]))
                self.tc2.append(ltc2)
                # print(d2)

            for j in range(self.J):
                column = data.readline().split()
                ltc3 = []
                for w in range(self.W):
                    ltc3.append(float(column[w]))
                self.tc3.append(ltc3)

            for i in range(self.I):
                column = data.readline().split()
                lp1 = []
                for j in range(self.J):
                    lp1.append(float(column[j]))
                self.p1.append(lp1)

            for j in range(self.J):
                column = data.readline().split()
                lp2 = []
                for k in range(self.K):
                    lp2.append(float(column[k]))
                self.p2.append(lp2)

            for j in range(self.J):
                column = data.readline().split()
                lp3 = []
                for w in range(self.W):
                    lp3.append(float(column[w]))
                self.p3.append(lp3)

            for j in range(self.J):
                column = data.readline().split()
                lP1 = []
                for l in range(self.L):
                    lP1.append(float(column[l]))
                self.P1.append(lP1)

            for k in range(self.K):
                column = data.readline().split()
                lP2 = []
                for l in range(self.L):
                    lP2.append(float(column[l]))
                self.P2.append(lP2)

            for w in range(self.W):
                column = data.readline().split()
                lP3 = []
                for l in range(self.L):
                    lP3.append(float(column[l]))
                self.P3.append(lP3)

            for j in range(self.J):
                column = data.readline().split()
                lc1 = []
                for l in range(self.L):
                    lc1.append(float(column[l]))
                self.c1.append(lc1)

            for k in range(self.K):
                column = data.readline().split()
                lc2 = []
                for l in range(self.L):
                    lc2.append(float(column[l]))
                self.c2.append(lc2)

            for w in range(self.W):
                column = data.readline().split()
                lc3 = []
                for l in range(self.L):
                    lc3.append(float(column[l]))
                self.c3.append(lc3)
            # print(self.c3)

            for n in range(self.N):
                llcn1 = []
                for j in range(self.J):
                    column = data.readline().split()
                    lcn1 = []
                    for l in range(self.L):
                        lcn1.append(float(column[l]))
                    llcn1.append(lcn1)
                self.cn1.append(llcn1)
            # print(cn1)

            for n in range(self.N):
                llcn2 = []
                for k in range(self.K):
                    column = data.readline().split()
                    lcn2 = []
                    for l in range(self.L):
                        lcn2.append(float(column[l]))
                    llcn2.append(lcn2)
                self.cn2.append(llcn2)

            for n in range(self.N):
                llcn3 = []
                for w in range(self.W):
                    column = data.readline().split()
                    lcn3 = []
                    for l in range(self.L):
                        lcn3.append(float(column[l]))
                    llcn3.append(lcn3)
                self.cn3.append(llcn3)

    def build(self):
        try:
            # define models for 'master' and 'sub'
            self.master = GRBPY.Model("master")
            self.sub = GRBPY.Model("sub")

            # disable log information
            self.master.setParam("OutputFlag", 0)
            self.sub.setParam("OutputFlag", 0)

            # use lazy constraints
            self.master.setParam("LazyConstraints", 1)

            # disable presolving in subproblem
            self.sub.setParam("Presolve", 0)

            # required to obtain farkas dual
            self.sub.setParam("InfUnbdInfo", 1)

            # use dual simplex
            self.sub.setParam("Method", 1)

            # construct master problem
            for j in range(self.J):
                for l in range(self.L):
                    self.y1[j, l] = self.master.addVar(lb=0.0, ub=1.0, vtype=GRBPY.GRB.BINARY)

            for k in range(self.K):
                for l in range(self.L):
                    self.y2[k, l] = self.master.addVar(lb=0.0, ub=1.0, vtype=GRBPY.GRB.BINARY)

            for w in range(self.W):
                for l in range(self.L):
                    self.y3[w, l] = self.master.addVar(lb=0.0, ub=1.0, vtype=GRBPY.GRB.BINARY)

            self.maxshipcost = self.master.addVar(0.0, GRBPY.GRB.INFINITY, 0.0, GRBPY.GRB.CONTINUOUS)

            for j in range(self.J):
                self.c11[j] = self.master.addConstr(sum(self.y1[j, l] for l in range(self.L)) <= 1)

            for k in range(self.K):
                self.c12[k] = self.master.addConstr(sum(self.y2[k, l] for l in range(self.L)) <= 1)

            for w in range(self.W):
                self.c13[w] = self.master.addConstr(sum(self.y3[w, l] for l in range(self.L)) <= 1)

            self.c01 = self.master.addConstr(sum(self.D[i] for i in range(self.I)) <= sum(
                self.C1[l] * self.y1[j, l] for j in range(self.J) for l in range(self.L)))

            self.c02 = self.master.addConstr(sum(self.D[i] for i in range(self.I)) <= sum(
                self.C2[l] * self.y2[k, l] for k in range(self.K) for l in range(self.L)) + \
                                             sum(self.C3[l] * self.y3[w, l] for w in range(self.W) for l in
                                                 range(self.L)))

            self.master.setObjective(
                sum(self.P1[j][l] * self.d1[l] * self.y1[j, l] for j in range(self.J) for l in range(self.L)) + \
                sum(self.P2[k][l] * self.d2[l] * self.y2[k, l] for k in range(self.K) for l in range(self.L)) + \
                sum(self.P3[w][l] * self.d3[l] * self.y3[w, l] for w in range(self.W) for l in range(self.L)) + \
                self.maxshipcost, GRBPY.GRB.MINIMIZE)

            # construct subproblem
            self.u10 = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)
            self.v10 = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for n in range(self.N):
                self.u1[n] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for n in range(self.N):
                self.v1[n] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for n in range(self.N):
                self.u2[n] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for n in range(self.N):
                self.v2[n] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for n in range(self.N):
                self.u3[n] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for n in range(self.N):
                self.v3[n] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for n in range(self.N):
                self.beta1[n] = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for n in range(self.N):
                self.beta2[n] = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for n in range(self.N):
                self.beta3[n] = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            self.beta4 = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)
            self.beta5 = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)
            self.beta6 = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            self.c = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for i in range(self.I):
                for j in range(self.J):
                    self.f1[i, j] = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for j in range(self.J):
                for k in range(self.K):
                    self.f2[j, k] = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for j in range(self.J):
                for w in range(self.W):
                    self.f3[j, w] = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for i in range(self.I):
                for j in range(self.J):
                    self.x1[i, j] = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for j in range(self.J):
                for k in range(self.K):
                    self.x2[j, k] = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for j in range(self.J):
                for w in range(self.W):
                    self.x3[j, w] = self.sub.addVar(lb=0.0, ub=GRBPY.GRB.INFINITY, vtype=GRBPY.GRB.CONTINUOUS)

            for i in range(self.I):
                self.alpha[i] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY,
                                                vtype=GRBPY.GRB.CONTINUOUS)

            for j in range(self.J):
                self.beta[j] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY,
                                               vtype=GRBPY.GRB.CONTINUOUS)

            for j in range(self.J):
                self.gamma[j] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=0, vtype=GRBPY.GRB.CONTINUOUS)

            for k in range(self.K):
                self.delt[k] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=0, vtype=GRBPY.GRB.CONTINUOUS)

            for w in range(self.W):
                self.epsilon[w] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=0, vtype=GRBPY.GRB.CONTINUOUS)

            for i in range(self.I):
                for j in range(self.J):
                    self.eta[i, j] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY,
                                                     vtype=GRBPY.GRB.CONTINUOUS)

            for j in range(self.J):
                for k in range(self.K):
                    self.theta[j, k] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY,
                                                       vtype=GRBPY.GRB.CONTINUOUS)

            for j in range(self.J):
                for w in range(self.W):
                    self.lamda[j, w] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=GRBPY.GRB.INFINITY,
                                                       vtype=GRBPY.GRB.CONTINUOUS)

            for j in range(self.J):
                for l in range(self.L):
                    self.Omega1[j, l] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=0, vtype=GRBPY.GRB.CONTINUOUS)

            for k in range(self.K):
                for l in range(self.L):
                    self.Omega2[k, l] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=0, vtype=GRBPY.GRB.CONTINUOUS)

            for w in range(self.W):
                for l in range(self.L):
                    self.Omega3[w, l] = self.sub.addVar(lb=-GRBPY.GRB.INFINITY, ub=0, vtype=GRBPY.GRB.CONTINUOUS)

            for i in range(self.I):
                self.c21[i] = self.sub.addConstr(sum(self.f1[i, j] for j in range(self.J)) == self.D[i])

            for j in range(self.J):
                self.c22[j] = self.sub.addConstr(sum(self.f1[i, j] for i in range(self.I)) - \
                                                 sum(self.f2[j, k] for k in range(self.K)) - sum(
                    self.f3[j, w] for w in range(self.W)) == 0)

            for j in range(self.J):
                self.c23[j] = self.sub.addConstr(
                    sum(self.f1[i, j] for i in range(self.I)) <= sum(self.C1[l] * 1 for l in range(self.L)))

            for k in range(self.K):
                self.c24[k] = self.sub.addConstr(
                    sum(self.f2[j, k] for j in range(self.J)) <= sum(self.C2[l] * 1 for l in range(self.L)))

            for w in range(self.W):
                self.c25[w] = self.sub.addConstr(
                    sum(self.f3[j, w] for j in range(self.J)) <= sum(self.C3[l] * 1 for l in range(self.L)))

            for i in range(self.I):
                for j in range(self.J):
                    self.c26[i, j] = self.sub.addConstr(self.f1[i, j] - self.Cbar1 * self.x1[i, j] == 0)

            for j in range(self.J):
                for k in range(self.K):
                    self.c27[j, k] = self.sub.addConstr(self.f2[j, k] - self.Cbar2 * self.x2[j, k] == 0)

            for j in range(self.J):
                for w in range(self.W):
                    self.c28[j, w] = self.sub.addConstr(self.f3[j, w] - self.Cbar2 * self.x3[j, w] == 0)

            for i in range(self.I):
                for j in range(self.J):
                    self.c31[i, j] = self.sub.addConstr(
                        self.alpha[i] + self.beta[j] + self.gamma[j] + self.eta[i, j] <= 0)

            for j in range(self.J):
                for k in range(self.K):
                    self.c32[j, k] = self.sub.addConstr(-self.beta[j] + self.delt[k] + self.theta[j, k] <= 0)

            for j in range(self.J):
                for w in range(self.W):
                    self.c33[j, w] = self.sub.addConstr(-self.beta[j] + self.epsilon[w] + self.lamda[j, w] <= 0)

            for i in range(self.I):
                for j in range(self.J):
                    self.c34[i, j] = self.sub.addConstr(-self.Cbar1 * self.eta[i, j] <= self.tc1[i][j])

            for j in range(self.J):
                for k in range(self.K):
                    self.c35[j, k] = self.sub.addConstr(-self.Cbar2 * self.theta[j, k] <= self.tc2[j][k])

            for j in range(self.J):
                for w in range(self.W):
                    self.c36[j, w] = self.sub.addConstr(-self.Cbar2 * self.lamda[j, w] <= self.tc3[j][w])

            for j in range(self.J):
                for l in range(self.L):
                    self.c411[j, l] = self.sub.addConstr(self.Omega1[j, l] - self.gamma[j] <= -self.M * 1 + self.M)

            for j in range(self.J):
                for l in range(self.L):
                    self.c412[j, l] = self.sub.addConstr(self.Omega1[j, l] - self.gamma[j] >= 0)

            for j in range(self.J):
                for l in range(self.L):
                    self.c413[j, l] = self.sub.addConstr(self.Omega1[j, l] >= -self.M * 1)

            for k in range(self.K):
                for l in range(self.L):
                    self.c421[k, l] = self.sub.addConstr(self.Omega2[k, l] - self.delt[k] <= - self.M * 1 + self.M)

            for k in range(self.K):
                for l in range(self.L):
                    self.c422[k, l] = self.sub.addConstr(self.Omega2[k, l] - self.delt[k] >= 0)

            for k in range(self.K):
                for l in range(self.L):
                    self.c423[k, l] = self.sub.addConstr(self.Omega2[k, l] >= -self.M * 1)

            for w in range(self.W):
                for l in range(self.L):
                    self.c431[w, l] = self.sub.addConstr(self.Omega3[w, l] - self.epsilon[w] <= - self.M * 1 + self.M)

            for w in range(self.W):
                for l in range(self.L):
                    self.c432[w, l] = self.sub.addConstr(self.Omega3[w, l] - self.epsilon[w] >= 0)

            for w in range(self.W):
                for l in range(self.L):
                    self.c433[w, l] = self.sub.addConstr(self.Omega3[w, l] >= -self.M * 1)

            self.c44 = self.sub.addConstr(sum(self.alpha[i] * self.D[i] for i in range(self.I)) + \
                                          sum(self.C1[l] * self.Omega1[j, l] for j in range(self.J) for l in
                                              range(self.L)) + \
                                          sum(self.C2[l] * self.Omega2[k, l] for k in range(self.K) for l in
                                              range(self.L)) + \
                                          sum(self.C3[l] * self.Omega3[w, l] for w in range(self.W) for l in
                                              range(self.L)) - \
                                          sum(self.tc1[i][j] * self.x1[i, j] for i in range(self.I) for j in
                                              range(self.J)) - \
                                          sum(self.tc2[j][k] * self.x2[j, k] for j in range(self.J) for k in
                                              range(self.K)) - \
                                          sum(self.tc3[j][w] * self.x3[j, w] for j in range(self.J) for w in
                                              range(self.W)) >= 0)

            self.c5 = self.sub.addConstr(
                self.u10 + self.c == sum(self.c1[j][l] * 1 for j in range(self.J) for l in range(self.L)) + \
                sum(self.c2[k][l] * 1 for k in range(self.K) for l in range(self.L)) + \
                sum(self.c3[w][l] * 1 for w in range(self.W) for l in range(self.L)))

            self.c6 = self.sub.addConstr(self.c <= self.C)

            for n in range(self.N):
                self.c51[n] = self.sub.addConstr(self.u1[n] + self.v1[n] ==
                                                 sum(self.cn1[n][j][l] * 1 for j in range(self.J) for l in
                                                     range(self.L)))

            for n in range(self.N):
                self.c52[n] = self.sub.addConstr(self.u2[n] + self.v2[n] ==
                                                 sum(self.cn2[n][k][l] * 1 for k in range(self.K) for l in
                                                     range(self.L)))

            for n in range(self.N):
                self.c53[n] = self.sub.addConstr(self.u3[n] + self.v3[n] ==
                                                 sum(self.cn3[n][w][l] * 1 for w in range(self.W) for l in
                                                     range(self.L)))

            self.c54 = self.sub.addConstr(self.u10 + sum(self.beta1[n] for n in range(self.N)) + \
                                          sum(self.beta2[n] for n in range(self.N)) + \
                                          sum(self.beta3[n] for n in range(self.N)) + \
                                          self.Gamma * (self.beta4 + self.beta5 + self.beta6) <= 0)

            for n in range(self.N):
                self.c551[n] = self.sub.addConstr(self.beta1[n] - self.u1[n] >= 0)

            for n in range(self.N):
                self.c552[n] = self.sub.addConstr(self.beta1[n] + self.u1[n] >= 0)

            for n in range(self.N):
                self.c561[n] = self.sub.addConstr(self.beta2[n] - self.u2[n] >= 0)

            for n in range(self.N):
                self.c562[n] = self.sub.addConstr(self.beta2[n] + self.u2[n] >= 0)

            for n in range(self.N):
                self.c571[n] = self.sub.addConstr(self.beta3[n] - self.u3[n] >= 0)

            for n in range(self.N):
                self.c572[n] = self.sub.addConstr(self.beta3[n] + self.u3[n] >= 0)

            for n in range(self.N):
                self.c581[n] = self.sub.addConstr(self.beta4 + self.sigma[n] * self.v1[n] >= 0)

            for n in range(self.N):
                self.c582[n] = self.sub.addConstr(self.beta4 - self.sigma[n] * self.v1[n] >= 0)

            for n in range(self.N):
                self.c591[n] = self.sub.addConstr(self.beta5 + self.sigma[n] * self.v2[n] >= 0)

            for n in range(self.N):
                self.c592[n] = self.sub.addConstr(self.beta5 - self.sigma[n] * self.v2[n] >= 0)

            for n in range(self.N):
                self.c5101[n] = self.sub.addConstr(self.beta6 + self.sigma[n] * self.v3[n] >= 0)

            for n in range(self.N):
                self.c5102[n] = self.sub.addConstr(self.beta6 - self.sigma[n] * self.v3[n] >= 0)

            self.sub.setObjective(
                sum(self.p1[i][j] * self.d * self.x1[i, j] for i in range(self.I) for j in range(self.J)) + \
                sum(self.p2[j][k] * self.d * self.x2[j, k] for j in range(self.J) for k in range(self.K)) + \
                sum(self.p3[j][w] * self.d * self.x3[j, w] for j in range(self.J) for w in range(self.W)),
                GRBPY.GRB.MINIMIZE)

        except GRBPY.GurobiError as e:
            print('Error code' + str(e.errno) + ': ' + str(e))
        except AttributeError as e:
            print('Encountered an attribute error: ' + str(e))

    def solve(self):
        # build 'master' and 'sub'
        self.build()

        # register callback
        self.master._iter = 0
        self.master._I = self.I
        self.master._J = self.J
        self.master._K = self.K
        self.master._W = self.W
        self.master._L = self.L
        self.master._N = self.N
        self.master._M = self.M
        self.master._C1 = self.C1
        self.master._C2 = self.C2
        self.master._C3 = self.C3
        self.master._Cbar1 = self.Cbar1
        self.master._Cbar2 = self.Cbar2
        self.master._D = self.D
        self.master._d = self.d
        self.master._d1 = self.d1
        self.master._d2 = self.d2
        self.master._d3 = self.d3
        self.master._tc1 = self.tc1
        self.master._tc2 = self.tc2
        self.master._tc3 = self.tc3
        self.master._p1 = self.p1
        self.master._p2 = self.p2
        self.master._p3 = self.p3
        self.master._P1 = self.P1
        self.master._P2 = self.P2
        self.master._P3 = self.P3
        self.master._c1 = self.c1
        self.master._c2 = self.c2
        self.master._c3 = self.c3
        self.master._cn1 = self.cn1
        self.master._cn2 = self.cn2
        self.master._cn3 = self.cn3
        self.master._c1 = self.c1
        self.master._c2 = self.c2
        self.master._c3 = self.c3
        self.master._sigma = self.sigma
        self.master._C = self.C
        self.master._Gamma = self.Gamma

        self.master._c = self.c
        self.master._u10 = self.u10
        self.master._v10 = self.v10
        self.master._u1 = self.u1
        self.master._v1 = self.v1
        self.master._u2 = self.u2
        self.master._v2 = self.v2
        self.master._u3 = self.u3
        self.master._v3 = self.v3
        self.master._beta1 = self.beta1
        self.master._beta2 = self.beta2
        self.master._beta3 = self.beta3
        self.master._beta4 = self.beta4
        self.master._beta5 = self.beta5
        self.master._beta6 = self.beta6

        self.master._c01 = self.c01
        self.master._c02 = self.c02

        self.master._c11 = self.c11
        self.master._c12 = self.c12
        self.master._c13 = self.c13
        self.master._c14 = self.c14

        self.master._c21 = self.c21
        self.master._c22 = self.c22
        self.master._c23 = self.c23
        self.master._c24 = self.c24
        self.master._c25 = self.c25
        self.master._c26 = self.c26
        self.master._c27 = self.c27
        self.master._c28 = self.c28

        self.master._c31 = self.c31
        self.master._c32 = self.c32
        self.master._c33 = self.c33
        self.master._c34 = self.c34
        self.master._c35 = self.c35
        self.master._c36 = self.c36

        self.master._c411 = self.c411
        self.master._c412 = self.c412
        self.master._c413 = self.c413
        self.master._c421 = self.c421
        self.master._c422 = self.c422
        self.master._c423 = self.c423
        self.master._c431 = self.c431
        self.master._c432 = self.c432
        self.master._c433 = self.c433
        self.master._c44 = self.c44

        self.master._c5 = self.c5
        self.master._c51 = self.c51
        self.master._c52 = self.c52
        self.master._c53 = self.c53
        self.master._c54 = self.c54
        self.master._c551 = self.c551
        self.master._c552 = self.c552
        self.master._c561 = self.c561
        self.master._c562 = self.c562
        self.master._c571 = self.c571
        self.master._c572 = self.c572
        self.master._c581 = self.c581
        self.master._c582 = self.c582
        self.master._c591 = self.c591
        self.master._c592 = self.c592
        self.master._c5101 = self.c5101
        self.master._c5102 = self.c5102

        self.master._c6 = self.c6

        self.master._y1 = self.y1
        self.master._y2 = self.y2
        self.master._y3 = self.y3
        self.master._maxshipcost = self.maxshipcost

        self.master._sub = self.sub

        # optimize master problem
        print("               *** Benders Decomposition Loop ***               ")
        self.master.optimize(cbwarehouse)
        print("                        *** End Loop ***                        ")

        self.sub.optimize()

    def report(self):
        print("               *** Summary Report ***               ")
        print("F1: %.6f" % (self.master.objval + self.sub.objval))
        print("F2: %.6f" % self.c.x)
        print("Variables:")
        print("               *** y ***               ")
        print([[self.y1[j,l].x for l in range(self.L)] for j in range(self.J)])
        print([[self.y2[k,l].x for l in range(self.L)] for k in range(self.K)])
        print([[self.y3[w,l].x for l in range(self.L)] for w in range(self.W)])
        print("               *** x ***               ")
        print([[self.x1[i, j].x for j in range(self.J)] for i in range(self.I)])
        print([[self.x2[j, k].x for k in range(self.K)] for j in range(self.J)])
        print([[self.x3[j, w].x for w in range(self.W)] for j in range(self.J)])
        print("               *** f ***               ")
        print([[self.f1[i, j].x for j in range(self.J)] for i in range(self.I)])
        print([[self.f2[j, k].x for k in range(self.K)] for j in range(self.J)])
        print([[self.f3[j, w].x for w in range(self.W)] for j in range(self.J)])

if __name__ == "__main__":
    warehouse = WareHouse()
    warehouse.read("d1")
    warehouse.solve()
    warehouse.report()
end = time.clock()
print('Running time: %s Seconds' % (end - start))