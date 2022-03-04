# FDM: a new package for numerical solution of ordinary differential equations in Sage

The Sage computer algebra system has very mediocre tools for the numerical integration of ordinary 
differential equations, but making computer experiments in it related to symbolic-numerical
calculations is very comfortable. The report presents a new package for numerical integration
of differential equations. We adhered to the following general principles: 
* numerical solutions are considered as elements of a new class, 
* the class definition provides tools for interpolating and visualizing the solution, 
* Richardson’s method for obtaining posterior estimates of errors is separate from the implementation of numerical methods, for which
two attributes – order of approximation and step or its analogue for quasi-equal grids are added in the class of numerical solutions. 

Our package focuses on visualizing the results of the calculation, including the construction of various kinds of auxiliary diagrams, including Richardson
diagrams. The implementation of Runge-Kutta method with arbitrary Butcher tableau,
for which a special class is introduced, is considered.

# Getting started
## Description of the initial problem
```
  sage: var("x1,x2,t")
  sage: problem1=InitialProblem([x2,-x1], [x1,x2], [0,1], 1)
```
## Description of the numerical solution
Explicit Runge-Kutta method with step dt=T/N:
```
  sage: P=erk(problem1, N=20)
```
Implicite Runge-Kutta method with step dt=T/N:
```
  sage: Q=irk(problem1, N=20)
```
Specifying a Butcher tableau:
```
  sage: B=Butcher_tableau(4,[[[1/4,1/4-1/6*sqrt(3)],[1/4+1/6*sqrt(3),1/4]],[1/2,1/2]], 'glrk','Gauss--Legendre RK methods')
  sage: R=irk(problem1, N=20,tableau=B)
```
## Iterpolation and plots
```
  sage: P.value(x1,pi) 
  sage: P.value(x1,pi,order=1)
  sage: P.plot(t,x1^2)
  sage: P.plot(x1,x2)
```
## Richardson estimate for error
```
  sage: L=[erk(problem1, N=20*2^n,tableau=butchers_list[1]) for n in range(10)]
  sage: richardson_plot(L,x1,9)
  sage: richardson(L[1],L[2],x1,9)
```
# References
* A. Baddour, M.D. Malykh, Richardson–Kalitkin method in abstract description, Discrete and Continuous Models and Applied Computational Science29 (3) (2021) 271–284.  DOI: [10.22363/2658-4670-2021-29-3-271-284](https://doi.org/10.22363/2658-4670-2021-29-3-271-284).
