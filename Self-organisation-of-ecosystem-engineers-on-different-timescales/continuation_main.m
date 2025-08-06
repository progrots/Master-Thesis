% Script for natural continuation of steady state solutions to
% Klausmeier_fast_pde. The script uses the function cont_Klausmeier.m, which uses the function findzeros.m, which
% uses the function newton.m.

clear variables
close all
clc

signtr = '12a-'; % Keeping track of different continuation experiments
%% Parameters

m = 1;
b = 0.1;
D = 0.01;
a = 3.5;
n = 1001;
L = 10*pi;
h = L/n;

%% Initial steady state + values of the fixed engineered construct (s0 + epsi * s)

% % 1 - spatially homogeneous initial conditions, from ode simulation
% sgn = append(signtr, 'ICh-');
% u0 = 1;
% v0 = 1;
% s0 = 1;
% t_range = [0,1000];
% options=odeset('RelTol',1e-10,'AbsTol',1e-12);
% [t1,ys] = ode15s(@(t,y) Klausmeier_plus_ode(t,y,a,m,b,1),t_range, [u0,v0,s0]);
% 
% u0 = ys(end,1); v0 = ys(end,2); s0st = ys(end,3);
% y0 = [u0 * ones(n,1); v0 * ones(n,1)]; s0 = s0st*ones(n,1);


% 2 - patterned initial conditions, from previous simulation
signtr = append(signtr, 'ICp-');
load("FullPDESimulation-12a-ac-ICh-\\y_end.mat","ys_end");
u0 = ys_end(1:n)'; v0 = ys_end(n+1:2*n)'; s = ys_end(2*n+1:end)';
y0 = ys_end(1:2*n)';

% 3 - manual
% signtr = append(signtr, 'ICm-');
% y0 = zeros(2*n,1);

% modification of initial conditions
% signtr = append(signtr, 'M-');
% y0 = [u0 * ones(n,1); v0 * ones(n,1)] + .1*cos((6*linspace(0,L,n)'*pi)/L)+.1

% 1 - parameter s from simulation; code works with engineered construct given by s0+epsi*s
signtr = append(signtr, 'ss-');
s0 = 0;
epsi = 1;
% s from simulation, assigned above
s0st = 0;
sst = s(1);

% 2 - parameter s as function; code works with soil state given by s0+epsi*s
% signtr = append(signtr, 'sf-');
% s_fun = @(x) cos((20*x*pi)/L+pi);
% sst = '@(x) cos((20*x*pi)/L+pi)';
% %s_fun = @(x) (alpha* cos((8*x*pi)/L+pi)+alpha).* ((x<2/8*L) + (x>4/8*L));%(alpha*cos(4*pi*x/L)+alpha); 
% %s_fun = @(x) (alpha*cos(6*pi*x/L)+alpha) .* ((x<1/6*L) + (x>3/6*L));
% s = s_fun(linspace(0,L,n)');
% % s0 from homogeneous steady state simulation above
% epsi = 1;

% % modification of s
% signtr = append(signtr, 'M-');
% % s(64:135) = zeros(136-64,1);
% s = s + 0.01*randn(n,1);


%% Newton iterations to find initial steady state
tol=10^(-12);
maxit=1000;
[y0s,iter]=newton(y0,@(y) Klausmeier_fast_pde(0, y, D, a, m, b, s0+epsi*s,n,h),@(y) JacKlausmeier_fast_pde(y, D, m, b, s0+epsi*s,n,h),tol,maxit,1);

%% Continuation
% cont_Klausmeier takes s0, epsi and s(x) and then continuates for soil state function s0 + epsi * s(x)
paramst = {'D', 'a', 'm', 'b', 'epsi','s0'};
param = [D,a,m,b,epsi,s0st]; % parameters needed for cont_Klausmeier (in that particular order)

% Continuation in parameter a, starting at steady state with vegetation
% and continuing down to a=0
p = 2; % parameter 'a'
signtr = sprintf('%sbifpar-%s',signtr, paramst{p});
par_start = a; % start in previously defined value for a
par_end = 0; % end at a=0
param(p) = par_start; % updating parameters
steps = 500; % number of continuation steps
[y_branch,eigvals,iters] = cont_Klausmeier(p, param, s, y0s, par_end, steps,h); % actual continuation
y_branch_singularities = y_branch(:,logical(eigvals(:,3))'); % retrieving solutions where the dominant eigenvalue changed sign


% Continuation in parameter a, starting at steady state with vegetation
% and continuing up to a=5
par_start2 = a;
par_end2 = 5;
param(p) = par_start2;
steps2 = 300;
[y_branch2,eigvals2,iters2] = cont_Klausmeier(p, param, s, y0s, par_end2, steps2,h);
y_branch_singularities2 = y_branch2(:,logical(eigvals2(:,3))');

% Continuation in parameter a, starting at barren state
par_start_null = 0;
par_end_null = 5;
param(p) = par_start_null;
steps0 = 800;
[y_branch_null,eigvals_null,iters_null] = cont_Klausmeier(p, param, s, zeros(2*n,1), par_end_null, steps0,h);
y_branch_singularities_null = y_branch_null(:,logical(eigvals_null(:,3))');

% Manual branch switching in case of branch point
% par_start_switched = eigvals(66,1);
% par_end_switched = 0;
% param(p) =  par_start_switched;
% y_init = y_branch(:,66);
% y_init(251:501) = y_init(251)*ones(251,1);
% y_init(751:end) = zeros(252,1);
% [y_branch_switched, eigvals_switched,iters_switched] = cont_Klausmeier(p,param,s,y_init,par_end_switched, steps2,h);


%% Plotting bifurcation diagrams
%Unpack variables
u_branch = y_branch(1:n,:);
u_branch2 = y_branch2(1:n,:);
u_branch_null = y_branch_null(1:n,:);
v_branch = y_branch(n+1:2*n,:);
v_branch2 = y_branch2(n+1:2*n,:);
v_branch_null = y_branch_null(n+1:2*n,:);

%Creating vectors of continuation parameter values
par2plot = linspace(par_start,par_end,steps+1);
par2plot2 = linspace(par_start2,par_end2,steps2+1);
par2plot0 = linspace(par_start_null,par_end_null,steps0+1);

%Splitting the positive solutions of v from the
%negative (nonsensical) solutions
index_pos_v = min(v_branch)>=0; 
index_neg_v = min(v_branch)<0;
v_pos = v_branch(:,index_pos_v);
v_neg = v_branch(:,index_neg_v);
u_pos = u_branch(:,index_pos_v);
u_neg = u_branch(:,index_neg_v);
par_pos_v = par2plot(index_pos_v);
par_neg_v = par2plot(index_neg_v);
%Retrieving bifurcation points
singu_pos_v = logical(eigvals(index_pos_v,3));
singu_neg_v = logical(eigvals(index_neg_v,3));

%Bifurcation diagram in the norm of u
bd_nrm_u = figure;
hold on

%Splitting stable from unstable solutions
stable_p = eigvals(index_pos_v,2)<=0;
unstable_p= eigvals(index_pos_v,2)>0;
stable_n = eigvals(index_neg_v,2)<=0;
unstable_n= eigvals(index_neg_v,2)>0;
%plotting first part of solutions
plot(par_pos_v(stable_p),sum(u_pos(:,stable_p))*h,'b') % plotting stable solutions in blue
plot(par_pos_v(unstable_p),sum(u_pos(:,unstable_p))*h,'r') % plotting unstable solutions in red
scatter(par_pos_v(singu_pos_v),sum(u_pos(:,singu_pos_v))*h,'r*') %marking bifurcation points with an asterisk
plot(par_neg_v,sum(u_neg)*h,'r--')
scatter(par_neg_v(singu_neg_v),sum(u_neg(:,singu_neg_v))*h,'r*')

%plotting second part of solutions
u_branch2 = y_branch2(1:n,:);
plot(par2plot2,sum(u_branch2)*h,'b');
scatter(par2plot2(logical(eigvals2(:,3))),sum(u_branch2(:,logical(eigvals2(:,3))))*h,'r*')

%Splitting stable from unstable solutions
u_branch_null = y_branch_null(1:n,:);
stable_null= eigvals_null(:,2)<=0;
unstable_null = eigvals_null(:,2)>0;
%plotting barren branch
plot(par2plot0(stable_null),sum(u_branch_null(:,stable_null))*h,'b');
plot(par2plot0(unstable_null),sum(u_branch_null(:,unstable_null))*h,'r');
scatter(par2plot0(logical(eigvals_null(:,3))),sum(u_branch_null(:,logical(eigvals_null(:,3))))*h,'r*')

% % In case of (manual) branch switching: plotting switched branch
% u_branch_switched = y_branch_switched(1:n,:);
% par2plot_switched = linspace(par_start_switched,0,steps2+1);
% stable_switched = eigvals_switched(:,2)<=0;
% unstable_switched = eigvals_switched(:,2)>0;
% plot(par2plot_switched(stable_switched),sum(u_branch_switched(:,stable_switched))*h,'b')
% scatter(par2plot_switched(unstable_switched),sum(u_branch_switched(:,unstable_switched))*h,'r.')
% singu_switched = logical(eigvals_switched(:,3));
% scatter(par2plot_switched(singu_switched),sum(u_branch_switched(:,singu_switched))*h, 'r*')

ylim([0,130]); ylabel('||u||_1')
xlim([0,3.5]); xlabel(paramst{p});

%Bifurcation diagram in the norm of v
bd_nrm_v = figure;
hold on

%plotting first part of solutions
plot(par_pos_v(stable_p),sum(v_pos(:,stable_p))*h,'b')
scatter(par_pos_v(unstable_p),sum(v_pos(:,unstable_p))*h,'r.')
scatter(par_pos_v(singu_pos_v),sum(v_pos(:,singu_pos_v))*h,'r*')
plot(par_neg_v,sum(v_neg)*h,'r--')
scatter(par_neg_v(singu_neg_v),sum(v_neg(:,singu_neg_v))*h,'r*')

%plotting second part of solutions
plot(par2plot2,sum(v_branch2)*h,'b');
scatter(par2plot2(logical(eigvals2(:,3))),sum(v_branch2(:,logical(eigvals2(:,3))))*h,'r*')

%plotting barren branch
plot(par2plot0(stable_null),sum(v_branch_null(:,stable_null))*h,'b');
scatter(par2plot0(unstable_null),sum(v_branch_null(:,unstable_null))*h,'r.');
scatter(par2plot0(logical(eigvals_null(:,3))),sum(v_branch_null(:,logical(eigvals_null(:,3))))*h,'r*')

% In case of (manual) branch switching: plotting switched branch
% v_branch_switched = y_branch_switched(n+1:end,:);
% plot(par2plot_switched(stable_switched),sum(v_branch_switched(:,stable_switched))*h,'b')
% scatter(par2plot_switched(unstable_switched),sum(v_branch_switched(:,unstable_switched))*h,'r.')
% scatter(par2plot_switched(singu_switched),sum(v_branch_switched(:,singu_switched))*h, 'r*')
% 
ylim([-0.01,90]); ylabel('||v||_1')
xlim([0,3.5]); xlabel(paramst{p});

%% Plotting solutions on branches

% Plotting steady state found by newton iterations from (modified) initial y0
stst = figure;
hold on
plot(y0s(1:n))
plot(y0s(n+1:2*n))
plot(s0+epsi*s)
ylim([0,2]); 
xlim([0,n]); xticks([0,n]); xticklabels({'0','10\pi'}); 
xlabel('x')
legend('u','v','s')

% spatial plots of u at several points on the vegetation branch
u_branch_fig = figure;
hold on
plot(u_pos(:,1:20:end))
plot(u_branch2(:,1:20:end))
xlim([0,n]); xticks([0,n]); xticklabels({'0','10\pi'}); 
xlabel('x'); ylabel('u');
ylim([0,2])

% spatial plots of v at several points on the vegetation branch
v_branch_fig = figure;
hold on
plot(v_pos(:,1:20:end))
plot(v_branch2(:,1:20:end))
xlim([0,n]); xticks([0,n]); xticklabels({'0','10\pi'}); 
xlabel('x'); ylabel('v')
ylim([0,2])

% in case of branch switching: spatial plots of v on the switched branch
% v_branch_switched_fig = figure;
% hold on
% plot(v_branch_switched(:,1:20:end))
% xlim([0,n]); xticks([0,n]); xticklabels({'0','10\pi'}); 
% xlabel('x'); ylabel('v')
% ylim([0,5])

% % spatial plots of v at singularities
% singu = figure;
% hold on
% plot(y_branch_singularities(n+1:2*n,:))
% plot(y_branch_singularities2(n+1:2*n,:))
% xlim([0,n]);
% xticks([0,n]); xticklabels({'0','10\pi'})

%% Saving

contname = append('ContinuationRR-',signtr);
mkdir(contname)

% save parameters and discretisation
pars_disc = {'m', m; 'b', b; 'D', D; 'a', a; '\epsilon', epsi; 's0', s0st; 's', sst; 'L',L; 'n',n; 'h', h};
writecell(pars_disc,sprintf('%s\\pars_and_disc.txt',contname),'Delimiter','tab')
save(sprintf('%s\\s.mat', contname),'s')
% save raw values
save(sprintf('%s\\y_branch.mat',contname),'y_branch')
save(sprintf('%s\\eigvals.mat',contname),'eigvals')

save(sprintf('%s\\y_branch2.mat',contname),'y_branch2')
save(sprintf('%s\\eigvals2.mat',contname),'eigvals2')

save(sprintf('%s\\y_branch_null.mat',contname),'y_branch_null')
save(sprintf('%s\\eigvals_null.mat',contname),'eigvals_null')

% save raw figures and pngs
formatnsave_fig(stst,0,sprintf('\\%s\\initial_ss_continuation',contname))
formatnsave_fig(u_branch_fig,1,sprintf('\\%s\\u_branch',contname))
formatnsave_fig(v_branch_fig,1,sprintf('\\%s\\v_branch',contname))
formatnsave_fig(bd_nrm_u,0,sprintf('\\%s\\bd_nrm_u',contname))
formatnsave_fig(bd_nrm_v,0,sprintf('\\%s\\bd_nrm_v',contname))