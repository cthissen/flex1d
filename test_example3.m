close all; clear all; clc

%% Broken beam, distributed load
% set material parameters
g = 9.81;      % m/s^2
rho_m = 3340;  % kg/m^3, density of material below beam. Determines bouyancy force. 


% set flexure parameters
Te = 25e3; % elastic thickness (meters)
E = 100e9; % young's modulous (100 GPA, Pa = kg/m^2)
v = 0.25;  % Poissons ratio
D = (E*Te^3)/(12*(1-v^2)); % rigidity parameter
lambda = (g*rho_m/(4*D))^(1/4); %
k = g*rho_m; % "elastic foundation" in lithosphere problems is the bouyancy force


% define loading
Tsed = 3e3;    % m, thickness of sedimentary load
rho_sed = 2700;% kg/m^3, density of material that infills deflection (e.g. sediments)


% build load vector
%... load extends from qxLeft to qxRight with constant magnitude q0
dx = 1e2; % meters
qxLeft  = -50e3;  % meters
qxRight  = 100e3; % meters
x = qxLeft:dx:qxRight; % load vector
q0 = rho_sed*Tsed*g;
qx = q0 + zeros(size(x)); % load, kg/m^2 ? 

% calculate deflection
xSol = qxLeft:(1e3):300e3; % vector of locations where we want to calculate deflection
deflection = flex1d(x,qx,xSol,lambda,k,'broken',qxLeft); % calculated deflection

% plotting
hFig = figure(1); clf
subplot(2,1,1)
qPlot = [0,0,qx,0,0];
xPlot = [xSol(1),x(1),x,x(end),xSol(end)];
plot(xPlot/1e3,qPlot/1e6);
ylabel('load (MPa)')
xlabel('x (km)')
title(sprintf('Rectangular load, broken plate'))

subplot(2,1,2)
plot(xSol/1e3,-deflection/1e3)
ylabel('Deflection (km)')
xlabel('x (km)')




PlotOpts = setdefaultplottingopts;
PlotOpts.figureSize = 'fullPage';
publishfigure(hFig,PlotOpts);
publishfigure(hFig.Children(1),PlotOpts);
publishfigure(hFig.Children(2),PlotOpts);
savefigure_cjt(hFig,'github_ex3','-png')
