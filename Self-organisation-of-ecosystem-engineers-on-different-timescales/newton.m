function  [x, iter] = newton(x,f,jacf,tol,maxit)
%DESCRIPTION
    % Newton iterations
%INPUT
    % x - initial guess
    % f - function to find zero for
    % jacf - jacobian of f
    % tol - allowed tolerance
    % maxit - maximum iterations
%OUTPUT
    % x - result of Newton iterations (if converged, f(x)=0)
    % iter - number of required iterations

    iter  = 0;
    nrmdx = Inf;
    while ((nrmdx > tol) & (iter < maxit))
        iter = iter + 1;
        dx = -jacf(x)\f(x);
        x = x + dx;
        nrmdx = norm(dx);
    end
    if ((iter == maxit) & (nrmdx > tol))
        fprintf('WARNING in newton.m, maximum iterations reached but not converged yet\n')
    end
end