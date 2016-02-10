%% Infinite beam, distributed load
close all; clear all; clc

% set material parameters
g = 9.81; % m/s^2
rho_m = 3340;  %kg/m^3, density of material below beam. Determines bouyancy force. 

% set flexure parameters
Te = 5e3; % elastic thickness (meters)
E = 100e9; % young's modulous (100 GPA, Pa = kg/m^2)
v = 0.25;  % Poissons ratio
D = (E*Te^3)/(12*(1-v^2));
lambda = (g*rho_m/(4*D))^(1/4);
k = g*rho_m; % "elastic foundation" in lithosphere problems is the bouyancy force


% define loading
Tsed = 3e3;   % m, thickness of sedimentary load
rho_sed = 2700;% kg/m^3, density of material that infills deflection (e.g. sediments)



% build load vector
%... load extends from qxLeft to qxRight with constant magnitude q0
dx = 1e2; % meters
qxLeft  = -50e3;  % meters
qxRight  = 100e3; % meters
x = qxLeft:dx:qxRight; % load vector
q0 = rho_sed*Tsed*g;
qx = q0 + zeros(size(x)); % load, kg/m^2 ? 

% vector of locations where we want to calculate deflection
xSol = -300e3:(1e3):300e3;

% analytic solution for rectangular load (Watts, eq 3.32 - 3.34)
y = zeros(size(xSol));
for iX = 1:numel(xSol)
    Q = qx(1)/(2*k);
    a = abs(qxLeft -  xSol(iX)); % distance to left limit of load
    b = abs(qxRight - xSol(iX)); % distance to right limit of load
    Da = exp(-lambda*a)*cos(lambda*a);
    Db = exp(-lambda*b)*cos(lambda*b);
    if xSol(iX) < qxLeft
        % points left of load (8a or 3.33)
        y(iX) = Q*(Da-Db);
    elseif xSol(iX) > qxRight
        % points right of load (9a 3.34)
        y(iX) = -Q*(Da-Db);
    else 
        % points under load (7a or 3.33)
        y(iX) = Q*(2-Da-Db);
    end
end

% superposition solution
ySuper = flex1d(x,qx,xSol,lambda,k,'infinite');

% calculate error
errSuperRectangle = (1/numel(ySuper))*sum((ySuper - y).^2);
fprintf('Mean-squared error: %e\n',errSuperRectangle);

% calculate airy isostatic deflection:
airyfnc = @(b1) (Tsed-b1)*(rho_sed)/(rho_m-rho_sed)-b1;
airyDef = fzero(airyfnc,0);

hFig = figure(1); clf
hAx(1) = subplot(3,1,1);
qPlot = [0,0,qx,0,0];
xPlot = [xSol(1),x(1),x,x(end),xSol(end)];
plot(xPlot/1e3,qPlot/1e6);
xlim([-300,300])
ylim([0,1.5*q0/1e6])
ylabel('load (MPa)')
title(sprintf('Rectangular load, infinite plate'))

hAx(2) = subplot(3,1,2);
xlim([-300,300])
hold on
plot(xSol/1e3,-ySuper/1e3,'-k','LineWidth',8)
plot(xSol/1e3,-y/1e3,'LineWidth',2)

airyPlotY = [0,0,airyDef,airyDef,0,0];
airyPlotX = [xSol(1),x(1),x(1),x(end),x(end),xSol(end)];
plot(airyPlotX/1e3,-airyPlotY/1e3,'--k')
ylabel('Deflection (km)')
box on
arrow([-100,-airyDef/1e3],[x(1)/1e3,-airyDef/1e3],'Length',8,'Width',2);
hAiryT = text(-200,-airyDef/1e3,'Airy Isostacy Depth');
hAiryT.Parent = hAx(2);
hAiryT.HorizontalAlignment = 'center';
publishfigure(hAiryT);

hLeg = legend('  Analytical','  Numerical');
hLeg.Location = 'southeast';
hLeg.Box = 'off';
hLeg.FontSize = 16;

hAx(3) = subplot(3,1,3);
xlim([-300,300])
hold on
plot(xSol/1e3,(y-ySuper)*1000)
ylabel('Deflection difference (mm)')
xlabel('x (km)')
box on

% % Make pretty figure
% PlotOpts = setdefaultplottingopts;
% PlotOpts.figureSize = 'fullPage';
% publishfigure(hFig,PlotOpts);
% publishfigure(hAx(1),PlotOpts);
% publishfigure(hAx(2),PlotOpts);
% publishfigure(hAx(3),PlotOpts);
% 
% savefigure_cjt(hFig,'github_benchmark1','-png')
