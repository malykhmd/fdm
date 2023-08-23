################
# FDM ver. 1.15 #
################



##################
# Cremona scheme #
# Malykh, 2022   #
##################

def cremona_hat(f,x,xx,ring=QQ):
    K=PolynomialRing(ring,x)
    M=K(f).monomials()
    C=K(f).coefficients()
    L=0
    S=[i==(i+j)/2 for [i,j] in zip(x,xx)]
    SS=[i==j for [i,j] in zip(x,xx)]
    for [c,m] in zip(C,M):
        if m.degree()==2:
            if m.is_univariate():
                a=SR(m.variable(0))
                L=L+c*a*a.subs(SS)
            else: 
                a=SR(m.variable(0))
                b=SR(m.variable(1))
                L=L+c*(a*b.subs(SS)+a.subs(SS)*b)/2
        else:
            L=L+c*(SR(m).subs(S))  
    return L

def cremona_step(problem, dt, ring=QQ):
    [F,x,x0,T]=problem.list()
    xx=list(var(['x'+str(i) for i in x]))
    eqs = [j-i-cremona_hat(f,x,xx, ring=ring)*dt for [i,j,f] in zip(x,xx,F)]
    ans = solve(eqs,xx)[0]
    return [i.subs(ans) for i in xx]

def cremona_scheme(problem, N=10, field=QQ):
    [f,x,x0,T]=problem.list()
    K=FractionField(PolynomialRing(field,x))
    t0=0
    ans=[[t0]+x0]
    dt=T/N
    X=[K(g) for g in cremona_step(problem,dt,ring=QQ)]
    for n in range(N):
        t0=t0+dt
        S={K(i):j for [i,j] in zip(x,x0)}
        x0=[field(i.subs(S)) for i in X]
        ans.append([t0]+x0)
    return Numsol(ans,[t]+x,dt,2,problem)
    
def cremona_det(field,x,M,s):
    K=PolynomialRing(field,x)
    A=matrix(len(M),len(M), lambda i,j: field(SR(M[j]).subs({xx:QQ(ss) for [xx,ss] in zip([t]+x,s.list()[i])})))
    return A

def cremona_scheme_general(problem,dt,N,ring=QQ):
    [F,x,x0,T]=problem.list()
    r=cremona_step(problem, dt, ring=ring)
    S=[i==j for [i,j] in zip(x,x0)]
    SS=[S]
    for n in range(N-1):
        S=[i==(j.subs(S)).full_simplify() for [i,j] in zip(x,r)]
        SS.append(S)
    return SS
