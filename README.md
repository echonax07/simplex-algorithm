# README

A python3 implementation of the two-phase simplex algorithm for solving linear programs, as outlined in "Introduction to Algorithms" (CLRS). This is primarily just for demonstrative purposes, for practical purposes I would recommend using the linear program solvers in numpy/scipy.

# How to use

Given a linear program in standard form, i.e.

<p align="center">
  <img src="https://i.imgur.com/ZQegHSh.png">
</p>

form lists `A, b, c` from the coefficients to give as input to the `simplex` function. Which will either return a bounded optimal solution (as a list of floats), print that there is no bounded optimal solution or print that there are no feasible solutions.

# About

The simplex algorithm finds an optimal solution to the linear program by iteratively checking vertices of a convex polygonal region formed by the inequalities defining the linear program. It can be shown mathematically that the optimal solution, if it exists, will be one of these vertices. This algorithm works in two phases.

## Phase 1

The original linear program (LP) is passed to `initialize_simplex` which attempts to solve a modified version of the original LP to find an initial feasible solution (that is, a vertex on the convex polygon formed by the constraints). If no such vertex is found, then there is no solutions to original LP. If an initial feasible solution is found, this solution along with a modified, but equivalent form of the original LP is returned.

## Phase 2

Given the initial feasible solution and equivalent LP returned from `initialize_simplex`, the `pivot` function chooses which neighboring vertex of the feasible region to consider in the next iteration of `simplex`. Each time choosing a neighboring vertex which increases the objective function (whose coefficients are stored in `c`). The `simplex` function stops when all neighboring vertices result in a decrease in the objective function.
