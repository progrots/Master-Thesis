function jac=JacKlausmeier_plus_pde(y, D, m, b, tau,n,h)
% Jacobian of Klausmeier_plus_pde
    % Linear part - approximation of spatial derivatives:
    A = spdiags([[1 -2 1].*ones(2*n,1);zeros(n,3)],-1:1,3*n,3*n);
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
    s = y(2*n+1:3*n);
    j11 = -1 - v.*s;
    j12 = -u.*s;
    j13 = -u.*v;
    j21 = v.*s.*(1-b*v);
    j22 = -m+u.*s-2*b*u.*v.*s;
    j23 = u.*v.*(1-b*v);
    j32 = ones(n,1)/tau;
    j33 = -j32;
    zr = zeros(n,1);
    NL = spdiags([[j21;j32;zr],[j11;j22;j33],[zr;j12;j23],[zr;zr;j13]],[-n,0,n,2*n],3*n,3*n);
    
    jac = A + NL;
end