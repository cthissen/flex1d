%% Infinite beam, triangular distributed load
close all; clear all; clc

% set material parameters
g = 9.81; % m/s^2
rho_m = 3340;  %kg/m^3, density of material below beam. Determines bouyancy force. 

% set flexure parameters
Te = 25e3; % elastic thickness (meters)
E = 100e9; % young's modulous (100 GPA, Pa = kg/m^2)
v = 0.25;  % Poissons ratio
D = (E*Te^3)/(12*(1-v^2));
lambda = (g*rho_m/(4*D))^(1/4);
k = g*rho_m; % "elastic foundation" in lithosphere problems is the bouyancy force


% define loading
Tsed = 10e3;   % m, thickness of sedimentary load
rho_sed = 2700;% kg/m^3, density of material that infills deflection (e.g. sediments)



% triangular load of sediments that increases in thickness from 0 to 10 km 
dx = 1e3;      % meters
A  = -125e3;   % meters
B  =  75e3;   % meters
l = B-A;
x  = A:dx:B; % load vector
q0 = g*rho_sed*Tsed;
qx = (q0/l)*x - q0*A/l; % load q(x)

% vector of locations where we want to calculate deflection
xSol = -300e3:(1e3):300e3;

% analytic solution for triangular
y = zeros(size(xSol));
for iX = 1:numel(xSol)
    Q = q0/(4*lambda*l*k);
    a = abs(A - xSol(iX)); % distance to left limit of load
    b = abs(B - xSol(iX)); % distance to right limit of load
    
    Ca = exp(-lambda*a)*(cos(lambda*a)-sin(lambda*a));
    Cb = exp(-lambda*b)*(cos(lambda*b)-sin(lambda*b));
    Da = exp(-lambda*a)*cos(lambda*a);
    Db = exp(-lambda*b)*cos(lambda*b);
    
    if xSol(iX) < A
        % points left of load (11a)
        y(iX) = Q*(Ca - Cb - 2*lambda*l*Db);
    elseif xSol(iX) > B
        % points right of load (12a)
        y(iX) = Q*(Ca - Cb + 2*lambda*l*Db);
    else 
        % points under load (10a)
        y(iX) = Q*(Ca - Cb - 2*lambda*l*Db + 4*lambda*a);
    
    
    end
end

% superposition solution
ySuper = flex1d(x,qx,xSol,lambda,k,'infinite');

% calculate error
errSuperRectangle = (1/numel(ySuper))*sum((ySuper - y).^2);
fprintf('Mean-squared error: %e (km)\n',errSuperRectangle/1e3);

% calculate airy isostatic deflection:
airyfnc = @(b1,axn) (axn-b1)*(rho_sed)/(rho_m-rho_sed)-b1;
for iA = 1:numel(x)
    axn = (Tsed/l)*x(iA) - Tsed*A/l; 
    airyZ = @(b1) airyfnc(b1,axn);
    airyDef(iA) = fzero(airyZ,0);
end
% fprintf('Airy deflection: %3.2f (km)\n',airyDef/1e3);

hFig = figure(1); clf
hAx(1) = subplot(3,1,1);
qPlot = [0,0,qx,0,0];
xPlot = [xSol(1),x(1),x,x(end),xSol(end)];
plot(xPlot/1e3,qPlot/1e6);
ylabel('load (MPa)')
xlabel('x (km)')
title(sprintf('Triangular load, infinite plate'))

hAx(2) = subplot(3,1,2);
hold on
plot(xSol/1e3,-ySuper/1e3,'-k','LineWidth',8)
plot(xSol/1e3,-y/1e3,'LineWidth',2)

airyPlotY = [0,airyDef,0,0];
airyPlotX = [xSol(1),x,x(end),xSol(end)];
plot(airyPlotX/1e3,-airyPlotY/1e3,'--k')
ylabel('Deflection (km)')
xlabel('x (km)')
box on
hAiryT = text(-275,-airyDef(end)/1e3,'Airy Isostacy');
hAiryT.Parent = hAx(2);
hAiryT.HorizontalAlignment = 'left';
publishfigure(hAiryT);
arrow([hAiryT.Extent(1)+1.1*hAiryT.Extent(3),-airyDef(end)/1e3],[x(150)/1e3,-airyDef(150)/1e3],'Length',8,'Width',1);

hLeg = legend('  Analytical','  Numerical');
hLeg.Location = 'southeast';
hLeg.Box = 'off';
hLeg.FontSize = 16;

hAx(3) = subplot(3,1,3);
hold on
plot(xSol/1e3,(y-ySuper)*1000)
ylabel('Deflection difference (mm)')
xlabel('x (km)')
box on

% PlotOpts = setdefaultplottingopts;
% PlotOpts.figureSize = 'fullPage';
% publishfigure(hFig,PlotOpts);
% publishfigure(hAx(1),PlotOpts);
% publishfigure(hAx(2),PlotOpts);
% publishfigure(hAx(3),PlotOpts);
% 
% savefigure_cjt(hFig,'github_benchmark2','-png')