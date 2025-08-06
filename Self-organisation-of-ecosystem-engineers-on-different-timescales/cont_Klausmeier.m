function  [y_branch,eigvals,iters] = cont_Klausmeier(par_to_cont, pars, s, y0, par_end, steps,h)
% FUNCTION NAME:
%   cont_Klausmeier
%
% DESCRIPTION:
%   Natural continuation of one branch of solutions to Klausmeier_fast_pde
%   in the parameter specified by par_to_cont, where the 'variable' s in
%   Klausmeier_fast_pde is given by s0 + epsi * s here.
%
% INPUT:
%   par_to_cont - parameter to continue: D=1,a=2,m=3,b=4,epsi=5,s0=6
%   pars - initial values parameters on a branch: D=pars(1); a=pars(2); m=pars(3); b=pars(4); epsi=pars(5); s0 = pars(6)
%   s - spatially variable part in s0 + epsi * s
%   y0 - initial steady state solution corresponding to par
%   par_end - the steady state solution is continued from pars(par_to_cont) up and until par_end
%   steps - number of steps taken until par_end
%   h - distance between gridpoints in discretisation
%
% OUTPUT:
%   y_branch - matrix with in column i the solution y on the branch in the i-th continuation step
%   eigvals - matrix containing [continuation parameter value, corresponding dominant eigenvalue of solution, whether eigenvalue changed sign (0,1)] 
%   iters - number of Newton iterations needed
%
% ASSUMPTIONS AND LIMITATIONS:
%   The initial steady state (pars,y) is a regular point (so not a
%   bifurcation point); might not stay on the branch for solutions too close to branching points.

par_init = pars(par_to_cont); % initial value of continuation parameter
n = floor(length(y0)/2); % number of gridpoints in discretisation
stepsize = (par_end-par_init)/steps; % continuation step size
cont_vals = par_init:stepsize:par_end; % continuation parameter values
y_branch = zeros(2*n,steps+1); % preallocation of output
iters = zeros(steps+1,1); % preallocation of output
y_branch(:,1) = y0;

eigvals = zeros(steps+1,3); % preallocation of output
jac = JacKlausmeier_fast_pde(y_branch(:,1),pars(1),pars(3),pars(4),pars(6) + pars(5)*s,n,h); % jacobian matrix corresponding to initial steady state
% Finding the dominant eigenvalue of the initial steady state:
wrn = warning('error', 'MATLAB:eigs:NotAllEigsConverged'); % raise an error when no eigenvalues are found
try
    dom_eig = eigs(jac,1,'largestreal',SubspaceDimension=300); % Find the largest real eigenvalue
catch % unless that method does not converge; in that case
    warning("Increase Krylov Subspace dimension and/or tolerance to use largestreal option. Using smallestabs instead.")
    dom_eig= eigs(jac,1,'smallestabs'); % Find the eigenvalue closest to zero. If the initial steady state is stable, this should have the same effect.
end
warning(wrn);
eigvals(1,:) = [cont_vals(1),dom_eig,0];

% Continuation:
for i=2:steps+1
    pars(par_to_cont) = cont_vals(i); % update the parameter values with the new value of the continuation parameter
    [y_branch(:,i),iters(i)] = findzeros(y_branch(:,i-1),pars,s,n,h); % find the nearest equilibrium, i.e., y such that Klausmeier_fast_pde is zero
    jac = JacKlausmeier_fast_pde(y_branch(:,i),pars(1),pars(3),pars(4),pars(6)+pars(5)*s,n,h); % calculate the jacobian corresponding to that equilibrium
    try
        dom_eig = max(real(eigs(jac,20,dom_eig))); % Track the dominant eigenvalue by looking for the 20 eigs closest to the previous one and taking the maximum
    catch 
        warning('matrix might contain inf or nan')
        dom_eig = inf;
    end
    if sign(real(dom_eig)) == sign(real(eigvals(i-1,2))) % Check whether the sign of the dominant eigenvalue has changed
        eigvals(i,:) = [cont_vals(i),dom_eig,0];
    else
        eigvals(i,:) =  [cont_vals(i),dom_eig,1];
    end

end
