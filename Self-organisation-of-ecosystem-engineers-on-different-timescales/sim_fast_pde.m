% Script for simulation of Klausmeier_fast_pde (discretisation of approximation of extended Klausmeier model)

clear variables
close all
clc

plt=1; % make plots?
timelapseplot = 0; % make timelapse plot?
simulation_signature = '12b-again'; % Keeping track of different simulations
%% parameters and discretisation

% Parameters
m = 1;
b = 0.1;
D = 0.01;
a0 = 3.0;
tau_a = 1000;

% 1 - a is constant
simulation_signature = append(simulation_signature, 'ac-');
a = @(t) a0 + 0*t; % a is constant but defined as a function for easier implementation
ast = 0;

% 2 - a is a function
% simulation_signature = append(simulation_signature, 'af-');
% a = @(t) max(0,a0 - t/tau_a); % linear decay at time scale tau_a
% ast = '@(t) max(0,a0 - t/tau_a)'

% Discretisation
n = 1001;
L = 10*pi;
h = L/n;

% Simulation time
t_end = 1000;
t_range = [0,t_end];

%% Initial conditions and parameter s
% % 1 - arbitrary
% simulation_signature = append(simulation_signature, 'ICarb-');
% u0 = 2*ones(n,1);
% v0 = cos(100*pi/L*linspace(0,L,n))'+1;
% y0 = [u0;v0];

% % 2 - spatially homogeneous, from ode simulation
% simulation_signature = append(simulation_signature, 'ICh-');
% u0 = 1;
% v0 = 1;
% s0 = 1;
% options=odeset('RelTol',1e-10,'AbsTol',1e-12);
% [t1,ys] = ode45(@(t,y) Klausmeier_plus(t,y,a0,m,b,1),t_range, [u0,v0,s0]);
% u0 = ys(end,1); v0 = ys(end,2); s0st = ys(end,3);
% y0 = [u0 * ones(n,1); v0 * ones(n,1)]; s0 = s0st*ones(n,1);

% 3 - patterned, from previous simulation
simulation_signature = append(simulation_signature, 'ICp-');
load("FullPDESimulation-12b-ICh-again-tau1000-/ics.mat",'y0');
u0 = y0(1:n); v0 = y0(n+1:2*n); s0 = y0(2*n+1:end);
y0 = [u0;v0];
% modification of initial conditions
% simulation_signature = append(simulation_signature, 'ICM-');
% y0 = [u0 * ones(n,1); v0 * ones(n,1)] + .1*cos((6*linspace(0,L,n)'*pi)/L)+.1

% % 1 - parameter s arbitrary
% s = ones(n,1);
% s0st = 'v0';

% 2 - parameter s from simulation
simulation_signature = append(simulation_signature, 'ss');
s = v0; % from simulation

% % 3 - parameter s as function
% simulation_signature = append(simulation_signature, 'sf');
% alpha = .2;
% s_fun = @(x) alpha* cos((2*x*pi)/L+pi)+alpha;
% s0st = '@(x) alpha* cos((2*x*pi)/L+pi)+alpha';
% s = s_fun(linspace(0,L,n)');

% % modification of s
% simulation_signature = append(simulation_signature, 'sM-');
% s(125:211) = zeros(212-125,1);
% s(1:246) = zeros(246,1);
% s(256:501) = zeros(502-256,1);

%% Numerical integration

[t,ys] = ode15s(@(t,y) Klausmeier_fast_pde(t, y, D, a(t), m, b,s,n,h),t_range,y0);
ys_end = ys(end,:);
nrm = norm(Klausmeier_fast_pde(0,ys_end',D,a(t(end)),m,b,s,n,h)) % display norm to check if we have reached a steady state


%% unpacking variables and plotting

%unpacking variables
u = ys(:,1:n);
v = ys(:,n+1:2*n);

%plotting spatial solutions u,v at t=0
ic = figure;
hold on
plot(y0(1:n))
plot(y0(n+1:end))
plot(s)
ylim([0,5])
xlim([0,n])
xticks([0,n])
xticklabels({'0','10\pi'})
xlabel('x')
legend('u','v','s')
title('Initial conditions')

%optional: plot spatial solutions u,v from simulation as time lapse
if timelapseplot
    ts = figure;
    xlim([0,n]);
    for i=1:10:length(t)
        clf()
        xlim([0,n]);
        ylim([0,5]);
        hold on
        plot(u(i,:))
        plot(v(i,:))
        plot(s)
        hold off
        pause(.1)
    end
end

%plot spatial solutions u,v at end of simulation
stst = figure;
hold on
plot(u(end,:))
plot(v(end,:))
plot(s)
ylim([0,5])
xlim([0,n])
xticks([0,n])
xticklabels({'0','10\pi'})
xlabel('x')
legend('u','v','s')
title('Steady State')


%% Finding the 'actual' steady state by Newton iterations
tol=10^(-12);
maxit=1000;
yi = ys(end,:)'; %initial guess: end of the simulation
[yos,iter]=newton(yi,@(y) Klausmeier_fast_pde(0, y, D, a(t(end)), m, b, s,n,h),@(y) JacKlausmeier_fast_pde(y, D, m, b, s,n,h),tol,maxit,1);

%unpacking variables
uos = yos(1:n);
vos = yos(1:2*n);

%plotting steady state from Newton iterations
stst_nwt = figure;
hold on
plot(u(end,:))
plot(v(end,:))
plot(s)
ylim([0,5])
xlim([0,n])
xticks([0,n])
xticklabels({'0','10\pi'})
xlabel('x')
legend('u','v','s')
title('Steady State from Newton iterations')

%% Saving

simname = append('Simulation-',simulation_signature);
mkdir(simname)

% save params and discretisation
pars_disc = {'m', m; 'b', b; 'D', D; 'a0', a0; 'tau_a', tau_a; 'a', ast; 's', s0st;'dt', dt; 't_end', t_range(end); 'L',L; 'n',n; 'h', h};
writecell(pars_disc,sprintf('%s\\pars_and_disc.txt',simname),'Delimiter','tab')
save(sprintf("%s\\ics.mat",simname),'y0')
% save raw values
save(sprintf('%s\\y_values.mat',simname),'ys')
save(sprintf('%s\\t_values.mat',simname),'t')
save(sprintf('%s\\y_end_values.mat',simname),'ys_end')
save(sprintf('%s\\s.mat',simname),'s')
% save raw figures and pngs
if plt
    formatnsave_fig(ic,0,sprintf('\\%s\\IC',simname))
    formatnsave_fig(stst,0,sprintf('\\%s\\steadystate',simname))
    formatnsave_fig(stst_nwt,0,sprintf('\\%s\\steadystate_newton',simname))
end