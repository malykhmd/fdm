################
# FDM ver. 2.0 #
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
# Тут внесены правки!!!
    def subs(self, u, abc, field=False):
        [f,x,x0,T]=self.list()
        if len(abc)==len(x):
            S=[i==j for [i,j] in zip(x,abc)]
        else:
            S=[t==abc[0]]+[i==j for [i,j] in zip(x,abc[1:])]
        if type(u)==type([]): 
            if field==False:
                ans=[uu.subs(S) for uu in u]
            else:
                ans=[field(uu.subs(S)) for uu in u]            
        else:
            if field==False:
                ans=u.subs(S)
            else:
                ans=field(u.subs(S))
        return ans
    def diff(self, u):
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
            u=self.diff(u)
            ans=ans + 1/factorial(i)*u*(tau-t)^i
        return ans
    def curvature(self):
        [f,x,x0,T]=self.list()
        a=[1] + f
        b=[0] + [self.diff(i) for i in f]
        k=sum([(a[i]*b[j]-a[j]*b[i])^2 for i in range(len(a)) for j in range(len(a)) if i<j])/sum([a_^2 for a_ in a])^(3/2)
        return k
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
        field=(P[-1][-1]).parent()
        while P[n][0] < t0:
            n=n+1
        if abs(P[n-1][0]- t0) < abs(P[n][0]- t0):
            n=n-1
        s=[i==j for [i,j] in zip(self.variables, P[n])]
        if t0==P[n][0]:
            ans= u.subs(s)
        else: 
            ans=self.problem.taylor(u,self.order+1).subs(s).subs(tau=t0)
        return field(ans)
    def zeros(self,u):
        P=self.points
        nums=[m for m in range(len(P)-1) if u.subs([i==j for [i,j] in \
        zip(self.variables, P[m])])* u.subs([i==j for [i,j] in zip(self.variables, P[m+1])])<0]
        polys=[self.problem.taylor(u,self.order+1).subs([i==j for [i,j] in zip(self.variables, P[n])]) for n in nums]
        ans=[find_root(f,P[n][0],P[n+1][0]) for [f,n] in zip(polys,nums) ]
        return ans
    def spline(self,u,t0):
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
    def plot(self, u, v, axes_labels='automatic', linestyle='solid', color='blue', points=False):
        S=[[x_==p_ for [x_,p_] in zip(self.variables,p)] for p in self.points]
        if axes_labels=='automatic':
            labels = ['$'+str(latex(u))+'$','$'+str(latex(v))+'$']
        else:
            labels = axes_labels
        if self.size()<51 or points:
            ans = point([[u.subs(s),v.subs(s)] for s in S], axes_labels=labels, color=color)
        else: 
            ans = line([[u.subs(s),v.subs(s)] for s in S], axes_labels=labels, linestyle=linestyle, color=color)
        return ans
    def plot3d(self, u, v, w, axes_labels='automatic', linestyle='solid', color='blue', points=False):
        S=[[x_==p_ for [x_,p_] in zip(self.variables,p)] for p in self.points]
        if axes_labels=='automatic':
            labels = ['$'+str(latex(u))+'$','$'+str(latex(v))+'$','$'+str(latex(w))+'$']
        else:
            labels = axes_labels
        if self.size()<51 or points:
            ans = point([[u.subs(s),v.subs(s),w.subs(s)] for s in S], axes_labels=labels, color=color)
        else: 
            ans = line([[u.subs(s),v.subs(s),w.subs(s)] for s in S], axes_labels=labels, linestyle=linestyle, color=color)
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

def richardson(P1,P2,u,t1,delta=0):
    r=P1.order+delta
    h1=P1.h
    h2=P2.h
    u1=P1.value(u,t1)
    u2=P2.value(u,t1)
    c=(u1-u2)/(h1^r-h2^r)
    if h1>h2:
        ans=[u2,c*h2^r]
    else:
        ans=[u1,c*h1^r]
    return ans      
    
def richardson_zeros(P1,P2,u,num=0,delta=0):
    r=P1.order+delta
    h1=P1.h
    h2=P2.h
    u1=P1.zeros(u)[num]
    u2=P2.zeros(u)[num]
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
    
def richardson_plot_zeros(P, u, num=0, nmin=0, nmax=oo):
    r=P[-1].order
    L=[[P_.h, abs(P_.zeros(u)[num] - P[-1].zeros(u)[num])/(1-(P[-1].h/P_.h)^r)] for P_ in P[:-1]]
    g1=list_plot_loglog(L, axes_labels=['$h$','$|E(Z('+str(latex(u))+'))|$'])
    L=[[log(a,10),log(b,10)] for [a,b] in L]
    L=L[nmin:min(nmax,len(L))]
    var("x")
    [a,b]=mnk(L)
    ll='$y='+str(latex(a.n(digits=3)*x+b.n(digits=3))) +'$ for the root $t='+str(latex(P[-1].zeros(u)[num]))+'$'
    g2=plot_loglog(10^b*x^a,(x,min([P_.h for P_ in P[:-1]]),max([P_.h for P_ in P[:-1]])), legend_label=ll)
    return  g1 + g2

#Нужно переделать 
def mnk(P): 
    vars=var('a,b') 
    s=sum([(a*P[n][0]+b-P[n][1])^2 for n in range(len(P))]) 
    eqs=[diff(s,a)==0, diff(s,b)==0] 
    S=solve(eqs,vars)[0] 
    return [RR(a.subs(S)),RR(b.subs(S))]

