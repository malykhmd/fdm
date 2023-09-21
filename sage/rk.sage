################
# FDM ver. 2.0 #
################

    
#######################
# Runge--Kutta method #
#######################    
    
def erk(problem, N=10, tableau=Butcher_tableau(4,[[[0,0,0,0],[1/2,0,0,0],[0,1/2,0,0],[0,0,1,0]], [1/6,1/3,1/3,1/6]], 'rk4','Standard rk method'), field=RR):
    [f,x,x0,T]=problem.list()
    t0=0
    ans=[[t0]+x0]
    dt=field(T/N)
    a=tableau.a(field=field)
    b=tableau.b(field=field)
    c=tableau.c(field=field)
    for n in range(N):
        k=[problem.subs(f,[t0+c[0]*dt]+x0, field=field)]
        for m in range(1,tableau.number_of_stages()):
            L=[t0+c[m]*dt] + [x0_ + sum([a_*k__ for [a_,k__] in zip(a[m],k_)])*dt \
               for [x0_,k_] in zip(x0,zip(*k))]
            k.append(problem.subs(f, L, field=field))
        t0=t0+dt
        x0=[x0_ + sum([b_*k__ for [b_,k__]  in zip(b,k_)])*dt \
            for [x0_,k_] in zip(x0,zip(*k))]
        ans.append([t0]+x0)
    return Numsol(ans,[t]+x,dt,tableau.order(),problem)


def erk_adaptive(problem, h=0.1, tableau=Butcher_tableau(4,[[[0,0,0,0],[1/2,0,0,0],[0,1/2,0,0],[0,0,1,0]], [1/6,1/3,1/3,1/6]], 'rk4','Standard rk method'), field=RR):
    [f,x,x0,T]=problem.list()
    if type(f)!=type([]):
        f=[f]
        x=[x]
        x0=[x0]
    ans=[[0]+x0+[0]]
    a=tableau.a(field=field)
    b=tableau.b(field=field)
    c=tableau.c(field=field)
    t0=0
    K=problem.curvature()
    while t0<T:
        L=[x_==RR(x0_) for [x_,x0_] in zip(x,x0)]
        dt=h*(K.subs(L))^(-2/5)
        k=[[f_c.subs(L) for f_c in f]]
        for m in range(1,tableau.number_of_stages()):
            L=[x_==x0_ + sum([a_*k__ for [a_,k__] in zip(a[m],k_)])*dt \
               for [x_,x0_,k_] in zip(x,x0,zip(*k))]
            k.append([f_.subs(L) for f_ in f])
        x0=[x0_ + sum([b_*k__ for [b_,k__]  in zip(b,k_)])*dt \
            for [x0_,k_] in zip(x0,zip(*k))]
        t0=t0+dt
        ans.append([t0]+x0+[dt])
    return Numsol(ans,[t]+x,h,tableau.order(),problem)

def irk(problem, N=10, eps=10^-10, M=10^2, tableau=Butcher_tableau(2,[[[1/2]],[1]], 'midpoint','midpoint methods'), field=RR):
    [f,x,x0,T]=problem.list()
    t0=0
    ans=[[t0]+x0]
    dt=T/N
    a=tableau.a(field=field)
    b=tableau.b(field=field)
    c=tableau.c(field=field)
    s=tableau.number_of_stages()
    while t0<T:
        k=[problem.subs(f,[t0]+x0, field=field) for i in range(s)]
        delta = oo
        i=0
        while delta>eps:
            kk = [problem.subs(f,[t0 +c[m]*dt] + [x0_ + sum([a_*k__ for [a_,k__] in zip(a[m],k_)])*dt for [x0_,k_] in zip(x0,zip(*k))], field=field) for m in range(s)]
            delta=(matrix(kk)-matrix(k)).norm()
            if i>M:
                print('error: the simple iteration method does not converge')
                break
            i=i+1
            k=kk
        t0=t0+dt
        x0=[x0_ + sum([b_*k__ for [b_,k__]  in zip(b,k_)])*dt \
            for [x0_,k_] in zip(x0,zip(*k))]
        ans.append([t0]+x0)
    return Numsol(ans,[t]+x,dt,tableau.order(),problem)
    
def irk_adaptive(problem, h=10^-1, eps=10^-10, M=10^2, tableau=Butcher_tableau(2,[[[1/2]],[1]], 'midpoint','midpoint methods'), field=RR, v= False):
    [f,x,x0,T]=problem.list()
    t0=0
    ans=[[t0]+x0]
    a=tableau.a(field=field)
    b=tableau.b(field=field)
    c=tableau.c(field=field)
    s=tableau.number_of_stages()
    jac=jacobian(f,x)
    while t0<T:
        if v == True:
            print('t='+latex(t0))
        dt=field(h/jac.subs([t==t0]+[i==j for [i,j] in zip(x,x0)]).norm())
        k=[problem.subs(f,[t0]+x0, field=field) for i in range(s)]
        delta = oo
        i=0
        while delta>eps:
            kk = [problem.subs(f,[t0 +c[m]*dt] + [x0_ + sum([a_*k__ for [a_,k__] in zip(a[m],k_)])*dt for [x0_,k_] in zip(x0,zip(*k))], field=field) for m in range(s)]
            delta=(matrix(kk)-matrix(k)).norm()
            if i>M:
                print('error: the simple iteration method does not converge')
                break
            i=i+1
            k=kk
        t0=t0+dt
        x0=[x0_ + sum([b_*k__ for [b_,k__]  in zip(b,k_)])*dt \
            for [x0_,k_] in zip(x0,zip(*k))]
        ans.append([t0]+x0)
    return Numsol(ans,[t]+x,dt,tableau.order(),problem)


