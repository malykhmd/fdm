def RosanesMatrices(f,g):
    E1=matrix(3,3,[[0,0,1],[0,0,0],[-1,0,0]])
    E2=matrix(3,3,[[0,0,0],[0,0,1],[0,-1,0]])
    var('u,v,w')
    S=[x==u/w, y==v/w]
    F=field[u,v,w](f.subs(S)*w^2)
    A=E1+dt/2*QuadraticForm(field[u,v,w](F)).matrix()  # В QuadraticForm матрица отличается от нашей тридиции в 2 раза, поэтому пополам. 
    G=field[u,v,w](g.subs(S)*w^2)
    B=E2+dt/2*QuadraticForm(field[u,v,w](G)).matrix()
    return [A,B]

def ProjectivePoint(eqs,vars,field):
    S=solve(eqs,vars)[0]
    A=ProjectiveSpace(SR, 2)([p.subs(S) for p in vars])
    A.clear_denominators()
    A=ProjectiveSpace(field[dt], 2)([a.expand() for a in list(A)])
    return A

def BileniarForm(A,p,pp):
    return (matrix(1,3,list(pp))*A*matrix(3,1,list(p)))[0][0]

def mod_ham(H,x):
    f=[-H.diff(x[1]), H.diff(x[0])]
    mh=H+dt/3*matrix([H.diff(xx) for xx in x])*(1-dt/2*jacobian(f,x))^-1*matrix(2,1,f)
    return SR(mh[0][0]).full_simplify()
# Здесь используется такой порядок: x=[p,q]. Формула взята у Суриса. 