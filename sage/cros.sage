##########################
# Cros                   #
# Badour Ali, 6.11.2017. #
##########################

def cros(problem, N=10, S=1, field=RR):
    [f,x,x0,T]=problem.list()
    if type(f)!=type([]):
        f=[f]
        x=[x]
        x0=[x0]
    dt=T/N/S
    M=len(x)
    K=(matrix.identity(M) - (1+i)/2*dt*jacobian(f,x))^(-1)*matrix(f).transpose()
    ans=[[0]+x0]
    for n in range(N*S):
        x0 = [field((x0[l]+dt*K[l][0]).subs([i==j for [i,j] in zip(x,x0)]).real()) for l in range(M)]
        if (n+1)/S in ZZ: 
            ans.append([(n+1)*dt]+x0) 
    return Numsol(ans,[t]+x,dt,2,problem)
    
def eff_order(problem, u, N=50):
    @parallel
    def foo(SS):
        ans= cros(problem, N=N, S=SS)
        return ans
    L=list(foo([1, 2, 2^2]))
    t1= [pts[0] for pts in L[0][1].list() if L[0][0][0][0]==1]
    y1= [problem.subs(u,pts) for pts in L[0][1].list() if L[0][0][0][0]==1]
    y2= [problem.subs(u,pts) for pts in L[1][1].list() if L[1][0][0][0]==2]
    y3= [problem.subs(u,pts) for pts in L[2][1].list() if L[2][0][0][0]==2^2]
    s=[[t1[n], RR(ln(abs((y2[n]-y1[n])/(y3[n]-y2[n])))/ln(2))] for n in range(1,N) if y3[n]!=y2[n] and y1[n]!=y2[n]]
    return line(s, axes_labels=["$t$","order"],  legend_label="r="+latex(s[-1][1])) + plot(2, (0,problem.T), color="red", linestyle="--")


