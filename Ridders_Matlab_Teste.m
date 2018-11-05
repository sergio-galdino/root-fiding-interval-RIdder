format long 
intvalinit('DisplayInfsup') 
tol =1e-14              
f = @(x)x^3-1 
[root,iter]=Ridder(f,0.6,1.5,tol) 
diam(root) 
f(root)  