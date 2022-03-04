################
# FDM ver. 1.1 #
################

#RR=RealField(200)

class InitialProblem:
    def __init__(self, f, x, x0, T):
        self.f = f
        self.x = x
        self.x0 = x0
        self.T = T
    def list(self):
        ans = [self.f, self.x, self.x0, self.T]
        if type(self.f)!=type([]):
            f=[self.f]
            x=[self.x]
            x0=[self.x0]
        return ans
    def subs(self, u, abc):
        if len(x1)==len(self.x):
            S=[i==j for [i,j] in zip(self.x,abc)]
        else:
            S=[t==abc[0]]+[i==j for [i,j] in zip(self.x,abc[1:])]
        if type(u)==type([]): 
            ans=[uu.subs(S) for uu in u]
        else:
            ans=u.subs(S)
        return ans
    def D(self, u):
        if type(u)==type([]): 
            ans=[sum([diff(uu,i)*j for [i,j] in zip(self.x,self.f)]) + diff(uu,t) for uu in u]
        else:
            ans=sum([diff(u,i)*j for [i,j] in zip(self.x,self.f)]) + diff(u,t)
        return ans
    def latex(self):
        print("\\left \\{ \\begin{aligned} &")
        print("".join([r'\frac{d}{dt}'+latex(xx)+'='+latex(ff) + r', \quad ' for [xx,ff] \
                       in zip(self.x[:len(x)-1],self.f[:len(x)-1])])\
                      + r' \frac{d}{dt}'+latex(self.x[-1])+'='+latex(self.f[-1])+', \\\\ &')
        print("".join([latex(xx)+'(0)='+latex(xx0)+r', \quad ' for [xx,xx0] in zip(self.x[:len(x)-1],self.x0[:len(x)-1])])\
                       + latex(self.x[-1])+'(0)='+latex(self.x0[-1]))
        print("\\end{aligned} \\right. ")

class Numsol:
    def __init__(self, points, variables,h,order):
        self.points = points
        self.variables = variables
        self.h = h
        self.order = order
    def size(self):
        return len(self.points)
    def list(self):
        return self.points
    def values(self,u,t0):
        P=self.points
        n=0
        while P[n][0] < t0:
            n=n+1
        L=[]
        for i in range(n-2,n+3):
            if i>=0 and i< self.size():
                L.append(P[i])
        S=[[v_==p__ for [v_,p__] in zip(self.variables,p_)] for p_ in L]
        return spline([[t.subs(s),u.subs(s)] for s in S])
    def value(self,u,t0,order=0):
        if order == 0: 
            ans = self.values(u,t0)(t0)
        elif order == 1 or order == 2:
            ans = self.values(u,t0).derivative(t0, order=order)
        return ans
    def plot(self,u,v):
        S=[[x_==p_ for [x_,p_] in zip(self.variables,p)] for p in self.points]
        labels = ['$'+str(latex(u))+'$','$'+str(latex(v))+'$']
        if self.size()<51:
            ans = point([[u.subs(s),v.subs(s)] for s in S], axes_labels=labels)
        else: 
            ans = line([[u.subs(s),v.subs(s)] for s in S], axes_labels=labels)
        return ans
    def plot_dt(self):
        t=list(zip(*self.list()))[0]
        P=[[t[i],t[i+1]-t[i]] for i in range(len(t)-1)]
        labels = ['$t$','$dt$']
        if self.size()<51:
            ans = point(P, axes_labels=labels)
        else: 
            ans = line(P, axes_labels=labels)
        return ans


#########################################
# Error estimate by Richardson-Kalitkin #
#########################################  

def richardson(P1,P2,u,t1,order=0, delta=0):
    r=P1.order-delta
    h1=P1.h
    h2=P2.h
    u1=P1.value(u,t1,order=order)
    u2=P2.value(u,t1,order=order)
    c=(u1-u2)/(h1^r-h2^r)
    if h1>h2:
        ans=[u2,c*h2^r]
    else:
        ans=[u1,c*h1^r]
    return ans      

def richardson_plot(P,u,t1,order=0, nmin=0, nmax=oo):
    r=P[-1].order
    L=[[P_.h, abs(P_.value(u,t1,order=order) - P[-1].value(u,t1,order=order))/(1-(P[-1].h/P_.h)^r)] for P_ in P[:-1]]
    g1=list_plot_loglog(L, axes_labels=['$h$','$|E('+str(latex(u))+')|$'])
    L=[[log(a,10),log(b,10)] for [a,b] in L]
    L=L[nmin:min(nmax,len(L))]
    var("x")
    [a,b]=mnk(L)
    ll='$y='+str(latex(a.n(digits=3)*x+b.n(digits=3))) +'$'
    g2=plot_loglog(10^b*x^a,(x,min([P_.h for P_ in P[:-1]]),max([P_.h for P_ in P[:-1]])), legend_label=ll)
    return  g1 + g2

#Нужно переделать 
def mnk(P): 
    vars=var('a,b') 
    s=sum([(a*P[n][0]+b-P[n][1])^2 for n in range(len(P))]) 
    eqs=[diff(s,a)==0, diff(s,b)==0] 
    S=solve(eqs,vars)[0] 
    return [RR(a.subs(S)),RR(b.subs(S))]

#######################
# Runge--Kutta method #
#######################

def latex_zero(a):
    if a==0:
        ans=latex('')
    else:
        ans=latex(a)
    return ans

class Butcher_tableau:
    def __init__(self, n, tableau, short_name=[], notes=[]):
        self.n = n
        self.tableau = tableau
        self.short_name = short_name
        self.notes = notes
    def order(self):
        return self.n
    def a(self,field=RR):
        return [[field(a_cc) for a_cc in a_c] for a_c in self.tableau[0]]
    def b(self,field=RR):
        return [field(b_c) for b_c in self.tableau[1]]
    def number_of_stages(self):
        return len(self.b())
    def c(self,field=RR):
        return [sum([field(a__) for a__ in a_]) for a_ in self.tableau[0]]
    def latex(self,field=SR):
        a=self.a(field=field)
        b=self.b(field=field)
        c=self.c(field=field)
        n = len(b)
        print('\\begin{array}{c|'+'c'*n+'}')
        for i in range(n):
            print(latex_zero(c[i]) + " & " + " & ".join([latex_zero(a[i][j]) for j in range(n)]) + "\\\\")
        print("\\hline")
        print(" & " + " & ".join([latex_zero(b[j]) for j in range(n)]))
        print('\\end{array}')
    def test(self):
        a=self.a()
        b=self.b()
        n = len(b)
        [print('error in a matrice, row no. '+ latex(i)) for i in a if len(i)!=n]
        print('ok')

load("butchers_list.sage")

def erk(problem, N=10, tableau=butchers_list[0]):
    [f,x,x0,T]=problem.list()
    t0=0
    ans=[[t0]+x0]
    dt=RR(T/N)
    a=tableau.a(field=RR)
    b=tableau.b(field=RR)
    c=tableau.c(field=RR)
    for n in range(N):
        k=[problem.subs(f,[t0+c[0]*dt]+x0)]
        for m in range(1,tableau.number_of_stages()):
            L=[t0+c[m]*dt] + [x0_ + sum([a_*k__ for [a_,k__] in zip(a[m],k_)])*dt \
               for [x0_,k_] in zip(x0,zip(*k))]
            k.append(problem.subs(f,L))
        t0=t0+dt
        x0=[x0_ + sum([b_*k__ for [b_,k__]  in zip(b,k_)])*dt \
            for [x0_,k_] in zip(x0,zip(*k))]
        ans.append([t0]+x0)
    return Numsol(ans,[t]+x,dt,tableau.order())

def curvature(f,x):
    a=[1] + [f_ for f_ in f]
    b=[0] + [sum([(diff(F,x_)*f_) for [x_,f_] in zip(x,f)]) for F in f]
    k=sum([(a[i]*b[j]-a[j]*b[i])^2 for i in range(len(a)) for j in range(len(a)) if i<j])/sum([a_^2 for a_ in a])^(3/2)
    return k

def erk_adoptive(problem, h=0.1, tableau=butchers_list[0]):
    [f,x,x0,T]=problem.list()
    if type(f)!=type([]):
        f=[f]
        x=[x]
        x0=[x0]
    ans=[[0]+x0+[0]]
    a=tableau.a(field=RR)
    b=tableau.b(field=RR)
    c=tableau.c(field=RR)
    t0=0
    K=curvature(f,x)
    while t0<T:
        L=[x_==RR(x0_) for [x_,x0_] in zip(x,x0)]
        dt=h*(K.subs(L))^(-2/5)
#        print(dt)
        k=[[f_c.subs(L) for f_c in f]]
        for m in range(1,tableau.number_of_stages()):
            L=[x_==x0_ + sum([a_*k__ for [a_,k__] in zip(a[m],k_)])*dt \
               for [x_,x0_,k_] in zip(x,x0,zip(*k))]
            k.append([f_.subs(L) for f_ in f])
        x0=[x0_ + sum([b_*k__ for [b_,k__]  in zip(b,k_)])*dt \
            for [x0_,k_] in zip(x0,zip(*k))]
        t0=t0+dt
        ans.append([t0]+x0+[dt])
    return Numsol(ans,[t]+x,h,tableau.order())

def irk(problem, N=10, eps=10^-10, M=10^2, tableau=butchers_list[2]):
    [f,x,x0,T]=problem.list()
    t0=0
    ans=[[t0]+x0]
    dt=RR(T/N)
    a=tableau.a(field=RR)
    b=tableau.b(field=RR)
    c=tableau.c(field=RR)
    s=tableau.number_of_stages()
    while t0<T:
        k=[problem.subs(f,[t0]+x0) for i in range(s)]
        delta = oo
        i=0
        while delta>eps:
            kk = [problem.subs(f,[t0 +c[m]*dt] + [x0_ + sum([a_*k__ for [a_,k__] in zip(a[m],k_)])*dt for [x0_,k_] in zip(x0,zip(*k))]) for m in range(s)]
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
    return Numsol(ans,[t]+x,dt,tableau.order())

################
# Adams method #
################

def adams(problem, N=10, r=2):
    [f,x,x0,T]=problem.list()
    if type(f)!=type([]):
        f=[f]
        x=[x]
        x0=[x0]
    ans=[[0]+x0]
    dt=RR(T/N)    
    D=lambda F: sum([(diff(F,x[i])*f[i]) for i in range(len(x))]) + diff(F,t)
    F = [f]
    g=f
    for i in range(r):
        g=[D(f_) for f_ in g]
        F.append(g)
    for n in range(N):
        L=[x_==x0_ for [x_,x0_] in zip(x,x0)] + [t==n*dt]
        x0=[x0[i] + sum([1/factorial(j+1)*F[j][i].subs(L)*dt^(j+1) for j in range(len(F))]) for i in range(len(f))]
        ans.append([(n+1)*dt]+x0)
    return Numsol(ans,[t]+x,dt,r+1)

def adams_adoptive(problem, h=10^-1, r=2):
    [f,x,x0,T]=problem.list()
    if type(f)!=type([]):
        f=[f]
        x=[x]
        x0=[x0]
    t0=0
    ans=[[t0]+x0]
    D=lambda G: sum([(diff(G,x[i])*f[i]) for i in range(len(x))]) + diff(G,t)
    F = [f]
    g=f
    for i in range(r+1):
        g=[D(f_) for f_ in g]
        F.append(g)
    while t0<T:
        L=[x_==RR(x0_) for [x_,x0_] in zip(x,x0)] + [t==RR(t0)]
        dt=h*(1/sqrt(sum([(1/factorial(r+1)*g_.subs(L))^2 for g_ in g])))^(1/(r+1))
        x0=[x0[i] + sum([1/factorial(j+1)*F[j][i].subs(L)*dt^(j+1) for j in range(len(F)-1)]) for i in range(len(f))]
        t0=t0+dt
        ans.append([t0]+x0)
    return Numsol(ans,[t]+x,h,r)

##################
#Midpoint method #
##################

#Solving of the eq. x=x0+f(x)
def simple_iteration_method(f,x,x0, eps=10^-20):
    x1=x0
    S=[i==j for [i,j] in zip(x,x0)]
    x2=[i+j.subs(S) for [i,j] in zip(x0,f)]
    n=0
    if abs(sum([(i-j)^2 for [i,j] in zip(x1,x2)]))>eps^2:
        while abs(sum([(i-j)^2 for [i,j] in zip(x1,x2)]))>eps^2:
            n=n+1
            x1=x2
            S=[i==j for [i,j] in zip(x,x1)]
            x2=[i+j.subs(S) for [i,j] in zip(x0,f)]
            if n>10^3:
                print('break')
                break
    return x2

def mpm_max_step_size(f,x,x0):
    M=len(f)
    S=[i==j for [i,j] in zip(x,x0)]
    J=matrix(RR, M, M, lambda n,m: diff(f[n], x[m]).subs(S))
    return RR(1/J.norm())
    
def mpm(problem, N=10, eps=10^-20):
    [f,x,x0,T]=problem.list()
    if type(f)!=type([]):
        f=[f]
        x=[x]
        x0=[x0]
    t0=0
    ans=[[0]+x0]
    dt=RR(T/N)
    while t0<T:
        if dt > mpm_max_step_size(f,x,x0):
            dt = mpm_max_step_size(f,x,x0)
        S=[i==(i+j)/2 for [i,j] in zip(x,x0)]
        F=[i.subs(S)*dt for i in f]
        x0=simple_iteration_method(F,x,x0, eps=eps)
        t0=t0+dt
        ans.append([t0]+x0)
    return Numsol(ans,[t]+x,dt,2)    	
