function  [y,iter]=findzeros(y0,pars,s,n,h)
% DESCRIPTION
%   Finding zeros of Klausmeier_fast_pde using Newton iterations
% INPUT:
%   y0 - initial guess for Newton iterations
%   pars - initial values parameters on a branch: D=pars(1); a=pars(2); m=pars(3); b=pars(4); epsi=pars(5); s0 = pars(6)
%   s - spatially variable part in s0 + epsi * s
%   n - number of discretisation gridpoints
%   h - distance between gridpoints in discretisation
%
% OUTPUT:
%   y - solution closest to y0 such that Klausmeier_fast_pde(0, y, D, a, m, b, s,n,h) = 0
%   iter - number of required Newton iterations

tol=10^(-6); % tolerance
maxit=1000; % maximum iterations
D=pars(1); a=pars(2); m=pars(3); b=pars(4); epsi=pars(5); s0 = pars(6);
[y,iter] = newton(y0,@(y) Klausmeier_fast_pde(0, y, D,a,m,b,s0+epsi*s,n,h),@(y) JacKlausmeier_fast_pde(y, D, m, b, s0+ epsi*s,n,h),tol,maxit,0);
end
