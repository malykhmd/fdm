# FDM: a new package for numerical solution of ordinary differential equations in Sage

[Sage](https://www.sagemath.org/) has very mediocre tools for the numerical integration of ordinary 
differential equations, but making computer experiments in it related to symbolic-numerical
calculations is very comfortable. Here we present a new package for numerical integration
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
  sage: problem1=Initial_problem([x1,x2], [x2,-x1], [0,1], 1)
```
## Description of the numerical solution
Explicit Runge-Kutta method with step dt=T/N:
```
  sage: P=erk(problem1, N=20)
```
Implicite Runge-Kutta method with quasistep h:
```
  sage: irk_adaptive(problem1, h=1)
```
Specifying a Butcher tableau:
```
  sage: B=symplectic_tableau(2)
  sage: B[0].latex(field=AA)
  sage: irk_adaptive(problem1, h=1, eps=10^-10, M=10^2, tableau=B[0])
```
## Interpolation and plots
```
  sage: P.value(x1,pi) 
  sage: P.plot(t,x1^2)
  sage: P.plot(x1,x2)
```
## Richardson estimate for error
```
  sage: L=[erk(problem1, N=20*2^n) for n in range(10)]
  sage: richardson_plot(L,x1,9)
  sage: richardson(L[1],L[2],x1,9)
```
# Authors 
The software was written by students and employees of RUDN since 2017:
* [Ali Baddour](https://orcid.org/0000-0001-8950-1781) (Syria)
* Ananina Luis Antonio Gonzalez (Ecuador)
* [Mikhail Malykh](https://orcid.org/0000-0001-6541-6603) (Russia)
* [Yu Ying](https://orcid.org/0000-0002-4105-2566) (China)
* Polina S. Chusovitina (Russia)
* Shiwei Wang (China)

# References
* Peter Stone. [Maple worksheets on the derivation of Runge-Kutta schemes](http://www.peterstone.name/Maplepgs/RKcoeff.html)
* Baddour A., Malykh M.D., Panin A.A., Sevastianov L.A. Numerical determination of the singularity order of a system of differential equations // Discrete and Continuous Models and Applied Computational Science. - 2020. - Vol. 28. - N. 1. - P. 17-34. doi: [10.22363/2658-4670-2020-28-1-17-34](https://doi.org/10.22363/2658-4670-2020-28-1-17-34)
* Baddour A., Malykh M.D. Richardson-Kalitkin method in abstract description // Discrete and Continuous Models and Applied Computational Science. - 2021. - Vol. 29. - N. 3. - P. 271-284. doi: [10.22363/2658-4670-2021-29-3-271-284](https://doi.org/10.22363/2658-4670-2021-29-3-271-284).
* Ying Y. The symbolic problems associated with Runge-Kutta methods and their solving in Sage // Discrete and Continuous Models and Applied Computational Science. - 2019. - Vol. 27. - N. 1. - P. 33-41. doi: [10.22363/2658-4670-2019-27-1-33-41](https://doi.org/10.22363/2658-4670-2019-27-1-33-41)
* Ying Y., Baddour A., Gerdt V.P., Malykh M.D., Sevastianov L.A. On the Quadratization of the Integrals for the Many-Body Problem. Mathematics 2021, 9, 3208. doi:[10.3390/math9243208](https://doi.org/10.3390/math9243208)
* M.D. Malykh, P.S. Chusovitina, Implementation of the Adams methodfor solving ordinary differential equations in the Sage computer algebrasystem, Discrete and Continuous Models and Applied Computational Sci-ence 31 (2) (2023) 164–173. DOI: [10.22363/2658-4670-2023-31-2-164-173](https://doi.org/10.22363/2658-4670-2023-31-2-164-173)
