import sympy as sp
import numpy as np
from scipy.optimize import root, minimize
from core import ConvexProblem

class KKTSolver:
    def __init__(self, problem: ConvexProblem):
        self.problem = problem
        self.logs = []

    def _build_numeric_system(self, conds):
        funcs = [sp.lambdify(self.problem.all_vars, cond, 'numpy') for cond in conds]
        def system(z):
            return [f(*z) for f in funcs]
        return system

    def _build_numeric_system_2(self, conds, used_vars):
        funcs = [sp.lambdify(used_vars, cond, 'numpy') for cond in conds]
        def system(z):
            return np.array([f(*z) for f in funcs])
        return system


    def solve_root(self, x0=None, verbose=True):
        conds = self.problem.get_kkt_conditions(full = True)
        system = self._build_numeric_system(conds)
        if x0 is None:
            x0 = np.ones(len(self.problem.all_vars))
        sol = root(system, x0)
        if verbose:
            print("Solver status:", sol.message)
        if not sol.success:
            return None
        return {str(var): val for var, val in zip(self.problem.all_vars, sol.x)}

    def solve_active_set(self, x0=None, tol=1e-5, max_iter=20, verbose=True):
        active_set = []

        for iteration in range(max_iter):
            if verbose:
                print(f"\n--- Iteration {iteration} ---")
                print(f"Active Set: {active_set}")

            conds = self.problem.get_kkt_conditions(active_set, full = False)
            used_vars = self.problem.used_vars
            system = self._build_numeric_system_2(conds, used_vars)
            x0 = x0 if x0 is not None else np.ones(len(used_vars))

            sol = root(system, x0)
            if not sol.success:
                print("Active set solver failed.")
                return None

            z = sol.x
            x_vars = z[:len(self.problem.vars)]
            lambda_vals = z[len(self.problem.vars):len(self.problem.vars)+len(self.problem.g)]
            x_dict = {str(var): val for var, val in zip(self.problem.all_vars, z)}

            # Check dual feasibility
            violated_lambdas = [i for i in active_set if lambda_vals[i] < -tol]
            if violated_lambdas:
                if verbose:
                    print(f"Removing constraint(s): {violated_lambdas}")
                active_set = [i for i in active_set if i not in violated_lambdas]
                x0 = z  # Update guess
                continue

            # Check for new active constraints
            inactive = [i for i in range(len(self.problem.g)) if i not in active_set]
            g_funcs = [sp.lambdify(self.problem.vars, self.problem.g[i], 'numpy') for i in inactive]
            g_vals = [f(*x_vars) for f in g_funcs]
            to_activate = [inactive[i] for i, val in enumerate(g_vals) if abs(val) < tol]

            if not to_activate:
                if verbose:
                    print("Converged with active set.")
                return x_dict

            if verbose:
                print(f"Activating constraint(s): {to_activate}")
            active_set.extend(to_activate)
            x0 = z  # warm-start

        print("Active set did not converge in time.")
        return None


    def solve_interior_point(self, x0=None, mu=1.0, tau=0.1, tol=1e-6, max_iter=10, verbose=True):
        # Barrier method for g(x) <= 0, keep equalities h(x)=0
        vars = self.problem.vars
        # numeric functions
        f_func = sp.lambdify(vars, self.problem.f, 'numpy')
        g_funcs = [sp.lambdify(vars, gi, 'numpy') for gi in self.problem.g]
        h_funcs = [sp.lambdify(vars, hi, 'numpy') for hi in self.problem.h]
        # starting point
        if x0 is None:
            x0 = np.ones(len(vars))
        e = 1e-4

        #    """
        def feas_obj(xk):
            # sum of squared violations
            penalty = 0.0
            for gi in g_funcs:
                gv = gi(*xk) + e
                if gv > 0:
                    penalty += gv**2
            for hi in h_funcs:
                penalty += hi(*xk)**2
            return penalty

        # Only run Phase I if current x0 is not strictly feasible
        need_phase1 = any(gi(*x0) + e >= 0 for gi in g_funcs)
        if need_phase1:
            if verbose:
                print("Running Phase I to find feasible start...")
            # unconstrained minimize penalty
            res_p1 = minimize(feas_obj, x0, method='SLSQP', options={'ftol': tol, 'disp': False})
            if not res_p1.success or feas_obj(res_p1.x) > tol:
                if verbose:
                    print("Unable to find strictly feasible start, infeasible problem or bad initialization.")
                return None
            x0 = res_p1.x
            if verbose:
                print(f"Phase I found feasible x0 = {x0}")
        #     """
        """
        def phase1_objective(z):
            # z = [x0_0, x0_1, ..., x0_{n-1}, s]
            return z[-1]  # minimize s

        # Build constraints for SLSQP:
        cons = []
        # (a) equality h_j(x)=0
        if h_funcs:
            cons.append({
               'type': 'eq',
               'fun': lambda z, funcs=h_funcs: np.array([f(*z[:-1]) for f in funcs])
            })
        # (b) g_i(x)+ε + s <= 0 → SLSQP uses g(x)>=0, so we flip:
        for gi in g_funcs:
            cons.append({
               'type': 'ineq',
               'fun': lambda z, gi=gi: -(gi(*z[:-1]) + e + z[-1])
            })

        # initial guess [old x0, s0]
        z0 = np.concatenate([x0, [1.0]])  # start s=1

        res1 = minimize(phase1_objective, z0,
                method='SLSQP',
                constraints=cons,
                bounds=[(None, None)]*len(x0) + [(0, None)],  # s>=0
                options={'ftol': tol, 'disp': False})

        if not res1.success or res1.x[-1] > tol:
            if verbose:
               print("Phase I failed to find strictly feasible x0 (s* = %.3e)." % res1.x[-1])
            return None

        # accept new x0
        x0 = res1.x[:-1]
        if verbose:
               print("Phase I found strictly interior x0 =", x0)


       """




        for i, gi in enumerate(g_funcs):
            val = gi(*x0)
            if val >= 0:
                if verbose:
                    print(f"Infeasible start for interior-point: g[{i}]= {val:.6f} >= 0")
                return None


        x = x0.copy()

        for iteration in range(max_iter):
            # build barrier objective
            def barrier_obj(xk):
                val = f_func(*xk)
                for gi in g_funcs:
                    gk = gi(*xk)
                    if gk >= 0:
                        return np.inf
                    val -= mu * np.log(-gk)
                return val
            # constraints: equality only
            cons = []
            if self.problem.h:
                cons = [{"type": "eq", "fun": (lambda xk, funcs=h_funcs: np.array([f(*xk) for f in funcs]))}]
            res = minimize(barrier_obj, x, method="SLSQP", constraints=cons, options={"ftol": tol, "disp": False})
            if not res.success:
                if verbose:
                    print(f"Barrier iter {iteration}: failed -", res.message)
                return None
            x = res.x
            if verbose:
                print(f"Barrier iter {iteration}: obj = {res.fun:.6f}, mu = {mu:.2e}")
            mu *= tau
        # final solution x only
        return {str(v): float(x[i]) for i, v in enumerate(vars)}

    def solve(self, method="root",**kwargs):
        if method == "root":
            return self.solve_root(**kwargs)
        elif method == "active_set":
            return self.solve_active_set(**kwargs)
        elif method == "interior_point":
            return self.solve_interior_point(**kwargs)
        else:
            raise ValueError("Unknown method. Use 'root', 'active_set', or 'interior_point'.")

