# ODE-Solver-Parallel-2
-This repository contains code that solves ODE with many boundary conditions that are auto generated based on boundary conditions entered by the user.  
-MUltifrontal Massively Parallel sparse direct Solver (MUMPS) has been used a solver for linear systems.  
-Message Passing interface (MPI) has been used while solving linear system each time.  
- **Only Solver part is in parallel and rest of the code runs in a serial fashion including matrix and right hand side(rhs) generation.**
