pkg load interval %  load interval package  
% Ridder    Finds  single  roots  of  a  function  in  given  range.
%
%           Ridder(func,a,b,tol)
%
%           Uses an interval version of Ridderís method to provide
%           rigorous  bounds  on  the single  roots  of  a  function  f.
%           Bounds are displayed as they are found.
%           Roots  are  displayed  if  radius  of  enclosure  <  tol
%           or if enclosure is no longer becoming tighter.   
%
%          INPUT:    func     function  of  nonlinear equation.
%                    x1,x2      limits of the interval containing the root.             
%                    tol      used as stopping criterion.
function [x,Iter] =Ridder(func,x1,x2,tol)
Iter=0;
x=infsup(x1,x2)
X=1.1*x;
f1=func(x1);
if f1 == 0; root = x1; return; end
f2=func(x2);
if f2 == 0; root = x2; return; end
if f1*f2 > 0
   'Root is not bracketed in (a,b)'
   x=infsup(NaN, NaN);
   return;
end
while or(wid(x) > tol , X < x)
    Iter=Iter+1;
    X=x;
    % Compute improved root from Ridderís formula
    x3 = mid(X); f3 = func(x3); 
    if f3 == 0; x = intval(x3); return; end
    s = sqrt(f3^2 - f1*f2);
    if s == 0; x=infsup(NaN, NaN); return; end
    dx = (x3 - x1)*f3/s;
    if (f1 - f2) < 0; dx = -dx; end
    x4 = x3 + dx; f4 = func(x4); 
    % Re-bracket the root 
    if f3*f4 > 0
       if f1*f4 < 0; x2 = x4; f2 = f4;
          else x1 = x4; f1 = f4;
       end
    else
        x1 = x3; x2 = x4; f1 = f3; f2 = f4; 
    end
    % Conversion to interval by infimum and supremum computed 
    % such that [x1,x2] is enclosed in interval x = infsup(x1,x2)
    if x1 < x2   
       x=infsup(x1, x2);
    else
       x=infsup(x2, x1);
    end 
 %  F=func(x); 
 %   if inf(F)*sup(F) > 0;x=X;return; end
end
endfunction