% Script for simulation of Klausmeier_plus_ode.

clear variables
close all
clc

simulation_signature = 'taus-bigger-taua-decaay'; % Keeping track of different simulations
plt=1; % make plots?
plot_nullclines = 0; % plot nullclines in the phase portrait?

%parameters
tau = 10000;
q = 0.5;
tau_a = tau^q;
dt = 0.01;
t_end = 500000;
% t_range = 0:dt:t_end;
a0 = 3;
a = @(t) max(a0 - t/tau_a,0*t); %a0 + 0*t; %a0+ exp(-t/2/tau_a).*a0.*cos(t/tau_a);%a0+exp(-t/tau_a); min(0 + t/tau_a,t./t);
ast = "@(t) max(a0 - t/tau_a,0*t);";
m = 1;
b = 0.1;
% when, a=0.5, m=0.1, h=1, u0=s0=0, v0=1, 
% vegetation dies out for tau = 100, whereas it persists for tau = 1

[t0,y0] = ode45(@(t,y) Klausmeier_plus_ode(t,y,a(0),m,b,tau),[0,4], [4,4,4]);

% u0 = 1;%1.2258;
% v0 = 1;%0.7947;
% s0 = 1;%0.7947;
y0 = y0(end,:);
% y0=[0.01,0.01,0.5]

options=odeset('RelTol',1e-7,'AbsTol',1e-5);
[t,y] = ode45(@(t,y) Klausmeier_plus_ode(t,y,a(t),m,b,tau),[0,t_end], y0,options);

%% Plotting
%defining nice plot colours
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
green = [78,180,0]/255;

%time series
ts = figure;
plot(t(2:end),y(2:end,:),'LineWidth', 1.2)
xlabel('t')
legend('u', 'v', 's')
ylim([0,4])

%phase plane
pp = figure;
hold on
plot3(y(:,1),y(:,2),y(:,3))
xlim([0,3]); ylim([0,3]); zlim([0,3])
xlabel('u')
ylabel('v')
zlabel('s')
if plot_nullclines
%nullclines
    [V,S] = meshgrid(linspace(0.05,5,100), linspace(0.05,5,100));
    U_udot = a0./(1+ V.*S);
    U_vdot = m./(S.*(1-b*V));
    s1=surf(U_udot,V,S, 'FaceAlpha',0.1);
    s2=surf(U_vdot,V,S, 'FaceAlpha',0.1);
    [U,V2] = meshgrid(linspace(0.05,5,100), linspace(0.05,5,100));
    s3=surf(U,V2,V2, 'FaceAlpha',0.1);
    s1.EdgeColor = 'none';
    s2.EdgeColor = 'none';
    s3.EdgeColor = 'none';
end

%plot simulation and slow manifolds: (v,s)-plot
slow2D = figure;
hold on
n=500;
S = linspace(0.05,1,n);
U = (a0*b + m)./(S+b);
V = 1/(a0/m*b+1)*(a0/m - 1./S);
plot(S,V,'b','LineWidth', 1.2)
plot(linspace(m/a0,1,n),zeros(n,1),'r','LineWidth', 1.2)
plot(linspace(0,m/a0,10),zeros(10,1),'b','LineWidth', 1)
plot(linspace(0,1,n),linspace(0,1,n),'k','LineStyle','--')
plot(y(:,3),y(:,2),'LineWidth', 2,color=[78,180,0]/255)
xlabel('s');ylabel('v');ylim([0,1]);xlim([0,1]);
set(gca,'xtick',[0,1])
set(gca,'ytick',[0,1])
%legend('M_v','M_0_-','M_0_+','v=s','simulation', Location='northwest')

%plot simulation and slow manifolds: (u,v,s)-plot
slow3d = figure;
hold on
plot3(U,V,S,'b','LineWidth', 1.2)
plot3(a0*ones(1,n),zeros(1,n),linspace(0,2,n),'r','LineWidth', 1.2)
plot3(y(:,1),y(:,2),y(:,3),color=[78,180,0]/255)
xlabel('u');ylabel('v');zlabel('s');ylim([0,1]);

%plot simulation in bifurcation diagram in a
bd = openfig('bd-a-v-1-01-1000.fig');
figure(bd)
hold on
plot(a(t),y(:,2), color=[78,180,0]/255,Linewidth=1.5)

%plot sim in fast system bifurcation diagram
fastbd = figure;
hold on
n=500;
s0=y0(3);
A = linspace(0.05,4,n);
V = 1./(A./m*b+1).*(A./m - 1./s0);
plot(A,V,'b','LineWidth', 1.2,'DisplayName','Stable manifold')
plot(linspace(m/s0,4,n),zeros(n,1),'r','LineWidth', 1.2, 'DisplayName','Unstable manifold')
plot(linspace(0,m/s0,10),zeros(10,1),'b','LineWidth', 1, 'DisplayName','Stable manifold')
% plot(linspace(0,4,n),linspace(0,4,n),'k','LineStyle','--')
plot(a(t),y(:,2),'LineWidth', 1.2,'Color',green, 'DisplayName','Simulation')
xlabel('a');ylabel('v');ylim([0,4]);
set(gca,'xtick',[0,4])
set(gca,'ytick',[0,4])
legend('','Unstable manifold', 'Stable manifold','Simulation')


%% Saving
simname = append('FullODESimulation-',simulation_signature);
mkdir(simname)

% save params, discretisation and initial condition)
pars_disc = {'m', m; 'b', b; 'a0', a0; 'tau_a', tau_a; 'tau',tau; 'a', ast; 'dt', dt; 't_end', t(end)};
writecell(pars_disc,sprintf('%s\\pars_and_disc.txt',simname),'Delimiter','tab')
save(sprintf("%s\\ics.mat",simname),'y0')
% save raw values
save(sprintf('%s\\y_values.mat',simname),'y')
save(sprintf('%s\\t_values.mat',simname),'t')
% save raw figures and pngs
if plt
    formatnsave_fig(pp,0,sprintf('\\%s\\phase-plane',simname))
    formatnsave_fig(bd,0,sprintf('\\%s\\sim-in-bd', simname))
    formatnsave_fig(fastbd,0,sprintf('\\%s\\sim-in-fast-bd', simname))
end