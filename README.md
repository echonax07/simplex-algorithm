# README

A python3 implementation of the two-phase simplex algorithm for solving linear programs, as outlined in "Introduction to Algorithms" (CLRS). This is primarily just for demonstrative purposes, for practical purposes I would recommend using the linear program solvers in numpy/scipy.

# How to use

Given a linear program in standard form, i.e.

![linprog](https://i.imgur.com/ZQegHSh.png)

form lists `A, b, c` from the coefficients to give as input to the simplex function. Which will either return a bounded optimal solution (as a list of floats), print that there is no bounded optimal solution or print that there are no feasible solutions.

# About

The simplex algorithm finds an optimal solution to the linear program by iteratively checking vertices of the convex polygonal region (feasible region) formed by the inequality. It can be shown mathematically that the optimal solution, if it exists, will be one of these vertices. This algorithm work into two phases.

## Phase 1

The original linear program (LP) is passed to `initialize_simplex` which uses the simplex algorithm on a modified version of the original LP to find an initial feasible solution (that is, a vertex on the convex polygon formed by the constraints). If no such vertex is found, there is no solutions to original LP.

## Phase 2

Given an initial feasible solution we iteratively move along neighboring vertices in the feasible region, each time choosing a neighboring vertex which increases the objective function (whose coefficients are stored in `c`).

A vertex is given by a solution of inequality constraints, it's neighboring vertices (i.e. vertices which share an edge) will share all but one of the same inequality constraints. The `pivot` function chooses an appropriate inequality (and hence a neighboring vertex) which increases the objective function.

The `simplex` function stops when no neighboring vertices can be found which result in an increase in the objective function.
