format long
tol =1e-14
f = @(x)(1+(1-5)^2)*x-(1-5*x)^2
[root,iter]=Ridder(f,0,0.2,tol)
wid(root)
f(root)