################
# FDM ver. 1.3 #
################

################
# Main classes #
################

class Initial_problem:
    def __init__(self, x, f, x0, T):
        self.f = f
        self.x = x
        self.x0 = x0
        self.T = T
    def list(self):
        if type(self.f)!=type([]):
            f=[self.f]
            x=[self.x]
            x0=[self.x0]
        else:
            f=self.f
            x=self.x
            x0=self.x0
        T=self.T
        return [f, x, x0, T]
    def subs(self, u, abc):
        [f,x,x0,T]=self.list()
        if len(abc)==len(x):
            S=[i==j for [i,j] in zip(x,abc)]
        else:
            S=[t==abc[0]]+[i==j for [i,j] in zip(x,abc[1:])]
        if type(u)==type([]): 
            ans=[uu.subs(S) for uu in u]
        else:
            ans=u.subs(S)
        return ans
    def D(self, u):
        [f,x,x0,T]=self.list()
        if type(u)==type([]): 
            ans=[sum([diff(uu,i)*j for [i,j] in zip(x,f)]) + diff(uu,t) for uu in u]
        else:
            ans=sum([diff(u,i)*j for [i,j] in zip(x,f)]) + diff(u,t)
        return ans
    def taylor(self,u,n):
        ans=u
        var('tau')
        for i in range(1,n+1):
            u=self.D(u)
            ans=ans + 1/factorial(i)*u*(tau-t)^i
        return ans
    def latex(self):
        [f,x,x0,T]=self.list()
        print("\\left \\{ \\begin{aligned} &")
        print("".join([r'\frac{d}{dt}'+latex(xx)+'='+latex(ff) + r', \quad ' for [xx,ff] \
                       in zip(x[:len(x)-1],f[:len(x)-1])])\
                      + r' \frac{d}{dt}'+latex(x[-1])+'='+latex(f[-1])+', \\\\ &')
        print("".join([latex(xx)+'(0)='+latex(xx0)+r', \quad ' for [xx,xx0] in zip(x[:len(x)-1],x0[:len(x)-1])])\
                       + latex(x[-1])+'(0)='+latex(x0[-1]))
        print("\\end{aligned} \\right. ")

class Numsol:
    def __init__(self, points, variables,h,order,problem):
        self.points = points
        self.variables = variables
        self.h = h
        self.order = order
        self.problem = problem
    def size(self):
        return len(self.points)
    def list(self):
        return self.points
    def value(self,u,t0):
        P=self.points
        n=0
        while P[n][0] < t0:
            n=n+1
        if abs(P[n-1][0]- t0) < abs(P[n][0]- t0):
            n=n-1
        s=[i==j for [i,j] in zip(self.variables, P[n])]
        return self.problem.taylor(u,self.order+1).subs(s).subs(tau=t0)
    def spline(self,u,t0):
        P=self.points
        n=0
        while P[n][0] < t0:
            n=n+1
        S=[[v_==p__ for [v_,p__] in zip(self.variables,p_)] for p_ in P[n-2:n+3]]
        return spline([[t.subs(s),u.subs(s)] for s in S])
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

def richardson_plot(P,u,t1, nmin=0, nmax=oo):
    r=P[-1].order
    L=[[P_.h, abs(P_.value(u,t1) - P[-1].value(u,t1))/(1-(P[-1].h/P_.h)^r)] for P_ in P[:-1]]
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
	
###################
# Butcher tableau #
###################

def latex_zero(a):
    if a==0:
        ans=latex('')
    elif a in AA or a in QQbar:
        ans=latex(a.radical_expression())
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
    def a(self,field=RR, radical_expression=False):
        if radical_expression==True:
            ans=[[field(a_cc).radical_expression() for a_cc in a_c] for a_c in self.tableau[0]]
        else: 
            ans=[[field(a_cc) for a_cc in a_c] for a_c in self.tableau[0]]
        return ans
    def b(self,field=RR, radical_expression=False):
        if radical_expression==True:
            ans=[field(b_c).radical_expression() for b_c in self.tableau[1]]
        else: 
            ans=[field(b_c) for b_c in self.tableau[1]]
        return ans
    def number_of_stages(self):
        return len(self.b())
    def c(self,field=RR, radical_expression=False):
        if radical_expression==True:
            ans=[sum([field(a__) for a__ in a_]).radical_expression() for a_ in self.tableau[0]]
        else: 
            ans=[sum([field(a__) for a__ in a_]) for a_ in self.tableau[0]]
        return ans
    def dic(self):
        s=self.number_of_stages()
        L=[self.tableau[0][i][j] for i in range(s) for j in range(s)]+self.tableau[1]
        S={i:L[i] for i in range(len(L))}
        return S
    def is_order(self,n):
        eqs=butcher_eqs(n,self.number_of_stages())
        L=[QQbar(eq.subs(self.dic())).minpoly() for eq in eqs]
        ans = prod([str(i)=='x' for i in L])==1
        return ans
    def is_symplectic(self):
        a=self.tableau[0]
        b=self.tableau[1]
        s=self.number_of_stages()
        L=[QQbar(b[i]*a[i][j]+b[j]*a[j][i]-b[i]*b[j]).minpoly() for i in range(s) for j in range(s)]
        ans = prod([str(i)=='x' for i in L])==1
        return ans
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

# From Yu Ying's phd these 
def butcher_list(p,s,implicit=True,symplectic=False):
    def lacroix(f, p):
        L=[f]
        D=lambda u: diff(u,x)*L[0]
        for  i in range(p-1):
            f=D(f)
            L.append(f)
        return L
    var('x,t,dt')
    b=var(['b'+str(n) for n in range(s)])
    A=matrix(s,s,var(['a'+str(n)+str(m) for n in range(s) for m in range(s)]))
    a=matrix(s,s,lambda i,j: (i>j or implicit)*A[i][j])
    c=[sum([a[i][j] for j in range(s)]) for i in range(s)]
    Kab=PolynomialRing(AA, list(a.variables())+list(b),order='lex')
    f=var(['f'+str(j) for j in range(p+1)])
    F=sum([f[j]*x^j/factorial(j) for j in range(p+1)])
    k=matrix(s,p+1,var(['k'+str(i)+str(j) for i in range(s) for  j in range(p+1)]))
    Slope=[sum([k[i][j]*dt^j/factorial(j) for j in range(p+1)]) for  i in range(s)]
    eqs=[Slope[i]-F.subs(x=dt*sum([a[i][j]*Slope[j] for j in range(s) if j<i or implicit])) for i in range(s)]
    sols=[]
    for j in range(p+1):
        sol=solve([eq.subs(dt=0)==0 for eq in eqs], [k[i][j] for i in range(s)])
        eqs=[diff(eq,dt).subs(sol) for eq in eqs]
        sols.append(sol)
    dx = dt*sum([b[i]*Slope[i].subs(sols) for i in range(s)]) 
    ser_dx=dx.series(dt, order=p+1).list()
    eqs = [lacroix(F, p)[j].subs(x=0) - factorial(j+1)*ser_dx[j+1] for  j in range(0,p) ]
    Kf=PolynomialRing(Kab, f)
    ans=[]
    for  j in range(0,p):
        ans+=Kf(eqs[j]).coefficients()
    if symplectic==True:
        ans=ans+[b[i]*a[i][j]+b[j]*a[j][i]-b[i]*b[j] for i in range(s) for j in range(s)]
# ans -- список уравнений, но из-за бага со словарем в variety его не получилось по другому скормить этому 
# методу, не получается передать переменные в subs. 
    V=(Kab*ans).variety(AA)
    T=[[[[Kab(a[i][j]).subs(v) for j in range(s)] for i in range(s)], [Kab(b[i]).subs(v) for i in range(s)]] \
       for v in V]
    return [Butcher_tableau(p,t) for t in T]
    
def butcher_eqs(p,s,implicit=True,symplectic=False):
    def lacroix(f, p):
        L=[f]
        D=lambda u: diff(u,x)*L[0]
        for  i in range(p-1):
            f=D(f)
            L.append(f)
        return L
    var('x,t,dt')
    b=var(['b'+str(n) for n in range(s)])
    A=matrix(s,s,var(['a'+str(n)+str(m) for n in range(s) for m in range(s)]))
    a=matrix(s,s,lambda i,j: (i>j or implicit)*A[i][j])
    c=[sum([a[i][j] for j in range(s)]) for i in range(s)]
    Kab=PolynomialRing(QQ, list(a.variables())+list(b),order='lex')
    f=var(['f'+str(j) for j in range(p+1)])
    F=sum([f[j]*x^j/factorial(j) for j in range(p+1)])
    k=matrix(s,p+1,var(['k'+str(i)+str(j) for i in range(s) for  j in range(p+1)]))
    Slope=[sum([k[i][j]*dt^j/factorial(j) for j in range(p+1)]) for  i in range(s)]
    eqs=[Slope[i]-F.subs(x=dt*sum([a[i][j]*Slope[j] for j in range(s) if j<i or implicit])) for i in range(s)]
    sols=[]
    for j in range(p+1):
        sol=solve([eq.subs(dt=0)==0 for eq in eqs], [k[i][j] for i in range(s)])
        eqs=[diff(eq,dt).subs(sol) for eq in eqs]
        sols.append(sol)
    dx = dt*sum([b[i]*Slope[i].subs(sols) for i in range(s)]) 
    ser_dx=dx.series(dt, order=p+1).list()
    eqs = [lacroix(F, p)[j].subs(x=0) - factorial(j+1)*ser_dx[j+1] for  j in range(0,p) ]
    Kf=PolynomialRing(Kab, f)
    ans=[]
    for  j in range(0,p):
        ans+=Kf(eqs[j]).coefficients()
    if symplectic==True:
        ans=ans+[b[i]*a[i][j]+b[j]*a[j][i]-b[i]*b[j] for i in range(s) for j in range(s)]
    return ans
    
#######################
# Runge--Kutta method #
#######################    
    
def erk(problem, N=10, tableau=Butcher_tableau(4,[[[0,0,0,0],[1/2,0,0,0],[0,1/2,0,0],[0,0,1,0]], [1/6,1/3,1/3,1/6]], 'rk4','Standard rk method')):
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
    return Numsol(ans,[t]+x,dt,tableau.order(),problem)

def curvature(f,x):
    a=[1] + [f_ for f_ in f]
    b=[0] + [sum([(diff(F,x_)*f_) for [x_,f_] in zip(x,f)]) for F in f]
    k=sum([(a[i]*b[j]-a[j]*b[i])^2 for i in range(len(a)) for j in range(len(a)) if i<j])/sum([a_^2 for a_ in a])^(3/2)
    return k


def erk_adaptive(problem, h=0.1, tableau=Butcher_tableau(4,[[[0,0,0,0],[1/2,0,0,0],[0,1/2,0,0],[0,0,1,0]], [1/6,1/3,1/3,1/6]], 'rk4','Standard rk method')):
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
    return Numsol(ans,[t]+x,h,tableau.order(),problem)

def irk(problem, N=10, eps=10^-10, M=10^2, tableau=Butcher_tableau(2,[[[1/2]],[1]], 'midpoint','midpoint methods')):
    [f,x,x0,T]=problem.list()
    t0=0
    ans=[[t0]+x0]
    dt=T/N
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
    return Numsol(ans,[t]+x,dt,tableau.order(),problem)
    
def irk_adaptive(problem, h=10^-1, eps=10^-10, M=10^2, tableau=Butcher_tableau(2,[[[1/2]],[1]], 'midpoint','midpoint methods')):
    [f,x,x0,T]=problem.list()
    t0=0
    ans=[[t0]+x0]
    a=tableau.a(field=RR)
    b=tableau.b(field=RR)
    c=tableau.c(field=RR)
    s=tableau.number_of_stages()
    jac=jacobian(f,x)
    while t0<T:
        dt=RR(h/jac.subs([t==t0]+[i==j for [i,j] in zip(x,x0)]).norm())
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
    return Numsol(ans,[t]+x,dt,tableau.order(),problem)

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
    return Numsol(ans,[t]+x,dt,r+1,problem)

def adams_adaptive(problem, h=10^-1, r=2):
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
    return Numsol(ans,[t]+x,h,r,problem)
    
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
    return Numsol(ans,[t]+x,dt,r+1,problem)

def adams_adaptive(problem, h=10^-1, r=2):
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
    return Numsol(ans,[t]+x,h,r+1,problem)
   	
