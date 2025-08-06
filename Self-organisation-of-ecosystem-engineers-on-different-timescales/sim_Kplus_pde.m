% Script for simulation of Klausmeier_plus_pde (discretisation of the extended Klausmeier model).

clear variables
close all
clc

plt = 1; % make plots?
timelapseplot = 0; % make timelapse plot?
simulation_signature = '19b-again-'; % Keeping track of different simulations
%% parameters and discretisation

% Parameters
m = 0.2;
b = 0.5;
D = 0.01;
tau = 1e9; 

% 1 - a is constant
% simulation_signature = append(simulation_signature, 'ac-');
% % For reference, the SN point and Turing point are given by
guess = 0.7; FindingTuring % aSN and aT
% a0 = 3.5;
% tau_a = 1;
% a = @(t) a0 + 0*t; % a is constant but formulated as a function for easier implementation
% ast = a0;

% 2 - a is a function
simulation_signature = append(simulation_signature, 'af-');
tau_a = 1e7;
a0 = 0.7;
a = @(t) max(0,a0 - t/tau_a); % linear decay at time scale tau_a
ast = '@(t) max(0,a0 - t/tau_a)';

% Discretisation
n = 501;
L = 10*pi;
h = L/n;

% Simulation time
t_end = 1e5;
t_range = [0,t_end];

%% Initial conditions

% 1 - arbitrary
% simulation_signature = append(simulation_signature, 'ICarb-');
% u0 = ones(n,1);
% v0 = -cos(20*pi/L*linspace(0,L,n))'+1;
% s0 = v0;
% y0 = [u0;v0;s0];

% % 2 - from homogeneous steady state
% simulation_signature = append(simulation_signature, 'ICh-');
% u0 = 1;
% v0 = 1;
% s0 = 1;
% options=odeset('RelTol',1e-10,'AbsTol',1e-12);
% % NOTE: to calculate the correct ss I might as well use tau=1, since tau
% % does not influence the location of equilibria.
% [t1,y0s] = ode45(@(t,y) Klausmeier_plus_ode(t,y,a(0),m,b,1),t_range, [u0;v0;s0]);
% y0 = [ones(n,1) * y0s(end,1);ones(n,1)*y0s(end,2);ones(n,1)*y0s(end,3)];

% 3 - patterend steady state from previous simulation
simulation_signature = append(simulation_signature, 'ICp-');
load("FullPDESimulation-12a-ac-ICh-/y_end.mat",'ys_end');
y0 = ys_end;

% modification to get the system out of equilibrium and hopefully onto a patterned state
% simulation_signature = append(simulation_signature, 'ICM-');
% y0=2*y0;
% y0(n+1+floor(n/2):2*n) = y0(n+1+floor(n/2):2*n)/2;
% y0(2*n+1+floor(n/2):end) = y0(2*n+1+floor(n/2):end)/2; 
% % a0=.7; % NOTE: this modified a0 pops up in the parameter txt-file
% % a = @(t) a0 + 0*t; % whereas a in the txt-file is the old value defined in ast
% % y0 = y0 + [zeros(n,1);0.01*sin(10*2*pi/L*linspace(0,L,n))';zeros(n,1)]; % periodic disturbance in v 
% % y0 = y0 + [zeros(n,1);0.01*randn(n,1);zeros(n,1)]; % random disturbance in v
% % y0 = y0 + 0.01*randn(1,3*n); % random disturbance in u,v,s
% % y0 = 10*y0;



%% Simulation
iend = 500; % max number of simulation extensions
options = odeset('RelTol',1e-9,'AbsTol',1e-11);
sol = ode15s(@(t,y) Klausmeier_plus_pde(t, y, D, a(t), m, b, tau,n,h),t_range,y0,options); % first simulation

% extend simulation until the simulation has (seemingly) reached steady state 
% or until the max number of extensions is reached
for i = 1:iend
    i
    nrm = norm(Klausmeier_plus_pde(0, sol.y(:,end), D, a(sol.x(end)), m, b, tau,n,h))
    if nrm < 1e-12
        break
    end

    solext = odextend(sol,[],(i+1)*t_end);
    clear sol;
    sol = solext;
    clear solext;
end

ys = sol.y'; % retrieving the simulated y-values
ys_end = ys(end,:); % separately storing y at the end of simulation
t = sol.x; % retrieving the simulation time steps

%% Plotting Simulation

%unpacking variables
u = ys(:,1:n);
v = ys(:,n+1:2*n);
s = ys(:,2*n+1:3*n);

%defining nice plot colours
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
green = [78,180,0]/255;

%optional: plot spatial solutions u,v,s from simulation as time lapse
if timelapseplot
    ts = figure;
    
    for i=1:10:length(t)
        clf()
        annotation('textbox', [.2 .5 .3 .3],'String',sprintf('t=%.0f',t(i)),'FitBoxToText','on');
        xlim([0,n]);
        ylim([0,1.5])
        hold on
        plot(u(i,:),'Color',blue)
        plot(v(i,:),'Color',red)
        plot(s(i,:),'Color',yellow) % I interchanged s and v to see if they truly overlap, and they do
        hold off
        pause(1)
    end
    close all
end

% optional: plot solutions u,v,s
if plt
    %plot spatial solutions u,v,s at t=0
    ic = figure;
    hold on
    plot(y0(1:n),'Color',blue,'DisplayName','u')
    plot(y0(n+1:2*n),'Color',red,'DisplayName','v')
    plot(y0(2*n+1:end),'Color',yellow,'DisplayName','s')
    ylim([0,4])
    xlim([0,n])
    xticks([0,n])
    xticklabels({'0','10\pi'})
    xlabel('x')
    legend
    %plot spatial solutions u,v,s at end of simulation
    sim_end = figure;
    hold on
    plot(u(end,:),'Color',blue,'DisplayName','u')
    plot(v(end,:),'Color',red,'DisplayName','v')
    plot(s(end,:),'Color',yellow','DisplayName','s')
    ylim([0,1])
    xlim([0,n])
    xticks([0,n])
    xticklabels({'0','10\pi'})
    xlabel('x')
    legend
    %plot ||u||_2, ||v||_2,||s||_2 as function of time
    nrmplt = figure;
    hold on
    plot(t,sum(u,2)*h,'DisplayName','||u||_1')
    plot(t,sum(v,2)*h,'DisplayName','||v||_1')
    plot(t,sum(s,2)*h,'DisplayName','||s||_1')
    xlim([0,t(end)])
    ylim([0,6])
    legend
    xlabel('t')

    % plot spatial solutions u,v,s for a specific time step
    % specificplt = figure;
    % ii = 91;
    % han=handle(annotation('textbox', [.2 .5 .3 .3],'String',sprintf('t=%.0f',t(ii)),'FitBoxToText','on','FontSize',16));
    % han.pinAtAffordance(1)
    % xlim([0,n]);
    % ylim([0,1.5])
    % hold on
    % plot(u(ii,:),'Color',blue,'DisplayName','u')
    % plot(v(ii,:),'Color',red, 'DisplayName','v')
    % plot(s(ii,:),'Color',yellow,'DisplayName','s') 
    % legend
    % xlabel('x'); xticks([0,n]),xticklabels({'0','10\pi'})

    %if a is a function: plot (||v||_2,a) and plot v(x,a) as a heatplot
    heatplt = figure;
    tiledlayout(2,1)
    ax1=nexttile;
    grid = linspace(0,L,n);
    surf(a(t),grid,v',EdgeColor="none")
    colorbar
    xlim([0,0.7]);ylim([0,L]);zlim([0,inf]);
    xlabel('a');ylabel('x');zlabel('v');
    yticks([0,L]); yticklabels({'0','10\pi'});
    ax2=nexttile;
    plot(a(t),sum(v,2)*h,'DisplayName','||v||_1')
    xlim([0,0.7]),ylim([0,inf]);
    xlabel('a');ylabel('||v||_1');
end
%% Finding and plotting patterned steady state -- Newton method does not work because jacobian badly scaled
% tol=10^(-6);
% maxit=1000;
% yi = ys(end,:)';
% [yos,iter]=newton(yi,@(y) Klausmeier_plus_pde(0, y, D, a(0), m, b, tau,n,h),@(y) JacKlausmeier_plus_pde(y, D, m, b, tau,n,h),tol,maxit,1);
% 
% uos = yos(1:n);
% vos = yos(1:2*n);
% sos = yos(2*n+1:3*n);
% 
% stst_nwt = figure;
% hold on
% plot(u(end,:))
% plot(v(end,:))
% plot(s(end,:))
% ylim([0,1])
% xlim([0,501])
% xticks([0,501])
% xticklabels({'0','2\pi'})
% xlabel('x')
% legend('u','v','s')

%% Saving

simname = append('FullPDESimulation-',simulation_signature);
mkdir(simname)

% save params, discretisation and initial condition)
pars_disc = {'m', m; 'b', b; 'D', D; 'a0', a0; 'tau_a', tau_a; 'tau',tau; 'a', ast; 't_end', t(end); 'L',L; 'n',n; 'h', h};
writecell(pars_disc,sprintf('%s\\pars_and_disc.txt',simname),'Delimiter','tab')
save(sprintf("%s\\ics.mat",simname),'y0')
% save raw values
save(sprintf('%s\\y_values.mat',simname),'ys')
save(sprintf('%s\\t_values.mat',simname),'t')
save(sprintf('%s\\y_end.mat',simname),'ys_end')
% save raw figures and pngs
if plt
    formatnsave_fig(ic,0,sprintf('\\%s\\IC',simname))
    formatnsave_fig(sim_end,0,sprintf('\\%s\\steadystate',simname))
    % formatnsave_fig(nrmplt,0,sprintf('\\%s\\normplot',simname))
    formatnsave_fig(heatplt,0,sprintf('\\%s\\heatplot',simname))
    % formatnsave_fig(specificplt,0,sprintf('\\%s\\specificplot',simname))
end