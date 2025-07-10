# 🧠 General-Purpose Constrained Optimization Solver using KKT Conditions

This project implements a flexible, symbolic-to-numeric optimization solver for constrained convex problems. It supports solving using Karush-Kuhn-Tucker (KKT) conditions with multiple backends: root solving, active-set, and interior-point (barrier) methods.

---

## 🚀 Features

- ✅ Symbolic formulation of objective and constraints using **SymPy**
- ⚙️ Modular architecture: define problems in `core.py`, solve via `solver.py`
- 🔁 Supports:
  - Root-based KKT solving
  - Active set method for inequality constraints
  - Interior-point method with automatic **Phase I** feasibility recovery
- 🔧 Automatic Lagrangian construction and KKT condition generation
- 🧪 Logs and prints solver status, feasibility, convergence, and variable values

---

## 📁 Project Structure

| File        | Description |
|-------------|-------------|
| `core.py`   | Defines the `ConvexProblem` class, symbolic expressions, KKT logic |
| `solver.py` | Implements `KKTSolver` with multiple methods (`root`, `active_set`, `interior_point`) |
| `example.py` or notebook | Demonstrates usage on example convex problems |

---

## 🔍 How It Works

You define an optimization problem symbolically:

```python
from sympy import symbols
from core import ConvexProblem

x1, x2 = symbols('x1 x2')
f = x1**2 + x2**2
g = [x1 - 0.1, x2 - 0.1]
h = [x1 + x2 - 1]

problem = ConvexProblem(f, [x1, x2], g_list=g, h_list=h)

## Then solve it using your preferred method:

from solver import KKTSolver
solver = KKTSolver(problem)

solution = solver.solve(method="interior_point", verbose=True)
```
---

## ✨ Available Solver Methods


Method |	Description
---------------------
`root`	 |  Solves full KKT system using scipy.optimize.root
`active_set` |	Iteratively adds/removes constraints based on multipliers
`interior_point` |	Barrier method for inequality constraints with feasibility recovery
