function jac=JacKlausmeier_fast_pde(y, D, m, b, s, n, h)
% Jacobian of Klausmeier_fast_pde
    % Linear part - approximation of spatial derivatives:
    A = spdiags([1 -2 1].*ones(2*n,1),-1:1,2*n,2*n);
    A(1,2) = 2;
    A(n,n-1) = 2;
    A(n,n+1) = 0;
    A(n+1,n) = 0;
    A(n+1,n+2) = 2;
    A(2*n,2*n-1) = 2;
    A(n+1:2*n,:) = A(n+1:2*n,:) * D;
    A = A/(h^2);
    % Nonlinear part:
    u = y(1:n);
    v = y(n+1:2*n);
    j11 = -1 - v.*s;
    j12 = -u.*s;
    j21 = v.*s.*(1-b*v);
    j22 = -m+u.*s-2*b*u.*v.*s;
    zr = zeros(n,1);
    NL = spdiags([[j21;zr],[j11;j22],[zr;j12]],[-n,0,n],2*n,2*n);
    
    jac = A + NL;
end