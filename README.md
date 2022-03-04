# FDM: a new package for numerical solution of ordinary differential equations in Sage

The Sage computer algebra system has very mediocre tools for the numerical integration of ordinary 
differential equations, but making computer experiments in it related to symbolic-numerical
calculations is very comfortable. The report presents a new package for numerical integration
of differential equations. When creating it, we adhered to the following general principles: 
* Numerical solutions are considered as elements of a new class, 
* the class definition provides tools for interpolating and visualizing the solution, 
* Richardson’s method for obtaining posterior estimates of errors is separate from the implementation of numerical methods, for which
two attributes – order of approximation and step or its analogue for quasi-equal grids are added in the class of numerical solutions. 

Our package focuses on visualizing the results of the calculation, including the construction of various kinds of auxiliary diagrams, including Richardson
diagrams. The implementation of Runge-Kutta method with arbitrary Butcher tables,
for which a special class is established, is considered.

# Getting started
## Description of the initial problem
```
  sage: var("x1,x2,t")
  sage: problem1=InitialProblem([x2,-x1], [x1,x2], [0,1], 1)
```
## Description of the numerical solution

```
  sage: P=erk(problem1, N=20,tableau=butchers_list[1])
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
