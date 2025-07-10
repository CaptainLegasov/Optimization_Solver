import sympy as sp
import numpy as np


class ConvexProblem:
    def __init__(self, f_expr, vars, g_list=None, h_list=None):
        self.vars = vars
        self.f = f_expr
        self.g = g_list or []
        self.h = h_list or []
        self.lambdas = sp.symbols(f"lambda0:{len(self.g)}")
        self.mus = sp.symbols(f"mu0:{len(self.h)}")
        self.all_vars = self.vars + list(self.lambdas) + list(self.mus)
        self.used_vars = []

    def get_kkt_conditions(self, active_set=None, full = True):
        L = self.f

        active_set = active_set if active_set is not None else range(len(self.g))

        if full:
            # Add all lambdas, regardless of active status
            for i in range(len(self.g)):
                L += self.lambdas[i] * self.g[i]
        else:
            for i in active_set:
                L += self.lambdas[i] * self.g[i]

        for j in range(len(self.h)):
            L += self.mus[j] * self.h[j]


        stationarity = [sp.diff(L, v) for v in self.vars]
        primal_eq = self.h
        if full:
        # Include all inequality constraints and slackness
           # active_g = [self.g[i] if i in active_set else 0 for i in range(len(self.g))]
            active_g = []
            slackness = [self.lambdas[i] * self.g[i] for i in active_set]
        else:
        # Use only active inequality constraints
            active_g = [self.g[i] for i in active_set]
            #slackness = [self.lambdas[i] * self.g[i] if i in active_set else self.lambdas[i]
             #            for i in range(len(self.g))
            inactive = [i for i in range(len(self.g)) if i not in active_set]
            slackness = [self.lambdas[i] for i in inactive]
            active_lambdas = [self.lambdas[i] for i in active_set]
            self.used_vars = self.vars + active_lambdas+ [self.lambdas[i] for i in inactive] + list(self.mus)


        if full:
            return stationarity + primal_eq + active_g + slackness

        else:
            return stationarity + primal_eq + active_g + slackness
