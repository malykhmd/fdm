################
# FDM ver. 1.15 #
################


################
# Adams method #
################

def adams(problem, N=10, r=2, field=RR):
    [f,x,x0,T]=problem.list()
    if type(f)!=type([]):
        f=[f]
        x=[x]
        x0=[x0]
    ans=[[0]+x0]
    dt=T/N
    F = [f]
    g=f
    for i in range(r):
        g=[problem.diff(f_) for f_ in g]
        F.append(g)
    for n in range(N):
        L=[x_==field(x0_) for [x_,x0_] in zip(x,x0)] + [t==n*dt]
        x0=[x0[i] + field(sum([1/factorial(j+1)*F[j][i].subs(L)*dt^(j+1) for j in range(len(F))])) for i in range(len(f))]
        ans.append([(n+1)*dt]+x0)
    return Numsol(ans,[t]+x,dt,r+1,problem)

def adams_adaptive(problem, h=10^-1, r=2, field=RR):
    [f,x,x0,T]=problem.list()
    if type(f)!=type([]):
        f=[f]
        x=[x]
        x0=[x0]
    t0=0
    ans=[[t0]+x0]
    F = [f]
    g=f
    for i in range(r+1):
        g=[problem.diff(f_) for f_ in g]
        F.append(g)
    while t0<T:
        L=[x_==field(x0_) for [x_,x0_] in zip(x,x0)] + [t==t0]
        dt=field(h*(1/sqrt(sum([(1/factorial(r+1)*g_.subs(L))^2 for g_ in g])))^(1/(r+1)))
        x0=[x0[i] + field(sum([1/factorial(j+1)*F[j][i].subs(L)*dt^(j+1) for j in range(len(F)-1)])) for i in range(len(f))]
        t0=t0+dt
        ans.append([t0]+x0)
    return Numsol(ans,[t]+x,h,r+1,problem)
    

