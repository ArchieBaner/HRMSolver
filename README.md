# HRMSolver
High Resolution Mass Solver using ILP strategy

Compilation:
---
gfortran -c ilp_solver.f90

gfortran hrms_solver.f90 ilp_solver.o -o hrms_solver

./hrms_solver

Features:
---
✔ Reads element data from elements.dat 📄

✔ Uses a custom ILP solver (Branch & Bound) for integer solutions 🔢

✔ Outputs all possible molecular formulas with deviations 📊

✔ Fast and efficient, no external solver needed 💨
