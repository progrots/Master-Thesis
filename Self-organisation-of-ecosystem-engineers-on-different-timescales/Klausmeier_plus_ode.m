function y_dot = Klausmeier_plus_ode(t, y, a, m, b, tau)
% Implementation of extended Klausmeier model without diffusion terms
% INPUT
%   t - time
%   y - 2-dim vector of variables u,v
%   a, m, b, tau - model parameters

    u = y(1);
    v = y(2);
    s = y(3);
    y_dot = zeros(3,1);
    y_dot(1) = a - u - u.*v.*s; %u_dot
    y_dot(2) = -m*v + u.*v.*s.*(1-b*v); %v_dot
    y_dot(3) = (v-s)/tau; %s_dot
