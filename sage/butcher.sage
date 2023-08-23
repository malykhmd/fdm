################
# FDM ver. 2.0 #
################


###################
# Butcher tableau #
###################

def latex_zero(a):
    if a==0:
        ans=latex('')
    elif a in AA or a in QQbar:
        ans=latex(QQbar(a).radical_expression())
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
        return len(self.tableau[1])
    def c(self,field=RR):
        return [sum([field(a__) for a__ in a_]) for a_ in self.tableau[0]]
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
    def latex(self,field=AA):
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
#    def test(self):
#        a=self.a()
#        b=self.b()
#        n = len(b)
#        [print('error in a matrice, row no. '+ latex(i)) for i in a if len(i)!=n]
#        print('ok')

# From Yu Ying's phd these 
def symplectic_tableau(s):
    def lacroix(f, p):
        L=[f]
        D=lambda u: diff(u,x)*L[0]
        for  i in range(p-1):
            f=D(f)
            L.append(f)
        return L
    var('x,t,dt')
    p=2*s
    b=var(['b'+str(n) for n in range(s)])
    a=matrix(s,s,var(['a'+str(n)+str(m) for n in range(s) for m in range(s)]))
    c=[sum([a[i][j] for j in range(s)]) for i in range(s)]
    Kab=PolynomialRing(AA, list(a.variables())+list(b),order='lex')
    f=var(['f'+str(j) for j in range(p+1)])
    F=sum([f[j]*x^j/factorial(j) for j in range(p+1)])
    k=matrix(s,p+1,var(['k'+str(i)+str(j) for i in range(s) for  j in range(p+1)]))
    Slope=[sum([k[i][j]*dt^j/factorial(j) for j in range(p+1)]) for  i in range(s)]
    eqs=[Slope[i]-F.subs(x=dt*sum([a[i][j]*Slope[j] for j in range(s)])) for i in range(s)]
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
    

