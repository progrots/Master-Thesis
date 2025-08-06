function y_dot = Klausmeier_plus_pde(t, y, D, a, m, b, tau,n,h)
%Discretisation of the extended Klausmeier pde
% on a one-dimensional grid with n equally spaced gridpoints at distance h.
% INPUT
%   t - time
%   y - vector of spatially discretised variables u,v
%   D, a, m, b, tau - model parameters
%   n - number of equidistant gridpoints
%   h - distance between gridpoints

    % Second order central difference approximation of spatial derivatives:
    A = spdiags([1 -2 1].*ones(2*n,1),-1:1,2*n,2*n);
    A(1,2) = 2;
    A(n,n-1) = 2;
    A(n,n+1) = 0;
    A(n+1,n) = 0;
    A(n+1,n+2) = 2;
    A(2*n,2*n-1) = 2;
    A(n+1:2*n,:) = A(n+1:2*n,:) * D;

    u = y(1:n);
    v = y(n+1:2*n);
    s = y(2*n+1:3*n);
    g = u.*v.*s;
    uv_dot = A*[u;v]/(h^2) + [a.*ones(n,1);zeros(n,1)] - [u;m*v] + [-g;g.*(1-b*v)];
    y_dot = [uv_dot; (v-s)/tau];

