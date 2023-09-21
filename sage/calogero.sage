##########################
# Calogero               #
# Shiwei Wang, 21.9.2023 #
##########################

def calogero_problem(ics,n,T=1,b=-1):
    ics=[QQ(icss) for icss in ics]
    q=var(['q'+ str(j) for j in range(n)])
    p=var(['p'+ str(j) for j in range(n)])
    H = sum([j^2 for j in p])/2 + sum([b/(j - k)^2 for j in q for k in q if j!=k])/2
    eqs=[diff(H,j) for j in p] + [-diff(H,j) for j in q]
    return Initial_problem(list(q+p),eqs,ics,T)
    
def calogero_q(ics,n,b=-1,t=1):
    var('q')
    ics=[QQ(icss) for icss in ics]
    t=QQ(t)
    Q_t0 = matrix(QQbar,n,n, lambda j,k: ics[j]*(j==k))
    L_t0 = matrix(QQbar,n,n, lambda j,k: ics[n+j]*(j==k) + (1-(j==k))/(ics[j]-ics[k] + (j==k)))
    F = (Q_t0 + sqrt(-b)*t*L_t0 - q).det()
    ans=AA[q](F).roots(multiplicities=False)
    if len(ans)<len(ics)/2:
        print('после момента столкновения')
        anss=ans
    else:
        order=[]
        for m in range(len(ics)/2):
            o=0
            for n in range(len(ics)/2):
                if ics[m]>ics[n]:
                    o=o+1
            order.append(o)
        anss=[]
        for m in range(len(order)):
            anss.append(ans[order[m]])
    return anss
    
def calogero_curve(ics,n,b=-1):
    var('q')
    ics=[QQ(icss) for icss in ics]
    Q_t0 = matrix(QQbar,n,n, lambda j,k: ics[j]*(j==k))
    L_t0 = matrix(QQbar,n,n, lambda j,k: ics[n+j]*(j==k) + (1-(j==k))/(ics[j]-ics[k] + (j==k)))
    F = (Q_t0 + sqrt(-b)*t*L_t0 - q).det()
    return F

def calogero_solution_crash(ics,n,b=-1):
    F=calogero_curve(ics,n,b=-1)
    eqs=[F, diff(F,q)]
    K=PolynomialRing(QQ,[q,t], order='lex')
    J=K*eqs
    G=J.groebner_basis()[-1]
    R=QQ[t](G).roots(AA,multiplicities=False)
    ans=min([r for r in R if r>0])
    return ans

