Here is computed the equilibrium temperature distribution on a rectangular 2D plate with Dirichlet boundary conditions on the square boundaries. 
One of the simplest examples of implementation of the Laplacian in numerical simulations.

There are 2 versions:
1- Time evolution:
It starts with a certain temperature distribution, and at every time t it computes the t+dt new temperature distribution. After some time the program converges.
The interruption can be made conditional on the convergence in many ways. Here for simplicity (and inexperience at that time,in 2020) was just asked to go on for a constant number of cycles.

2- Linear system of n x m equations and unknown. Matrix A of size (n x m ; n x m +1 ) is gauss reduced
