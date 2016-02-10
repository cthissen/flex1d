%% Broken beam, distributed rectangular load
close all; clear all; clc

% set material parameters
g = 9.81; % m/s^2
rho_m = 3340;  %kg/m^3, density of material below beam. Determines bouyancy force. 

% set flexure parameters
Te = 20e3; % elastic thickness (meters)
E = 100e9; % young's modulous (100 GPA, Pa = kg/m^2)
v = 0.25;  % Poissons ratio
D = (E*Te^3)/(12*(1-v^2));
lambda = (g*rho_m/(4*D))^(1/4);
k = g*rho_m; % "elastic foundation" in lithosphere problems is the bouyancy force
plateBreakPoint = -50e3;

% define loading
Tsed = 3e3;   % m, thickness of sedimentary load
rho_sed = 2700;% kg/m^3, density of material that infills deflection (e.g. sediments)


% build load vector
%... load extends from qxLeft to qxRight with constant magnitude q0
dx = 1e2; % meters
qxLeft  = -50e3;  % meters
qxRight  = 100e3; % meters
l = qxRight - qxLeft; % total distance of load
x = qxLeft:dx:qxRight; % load vector
q0 = rho_sed*Tsed*g;
qx = q0 + zeros(size(x)); % load, kg/m^2 ? 

% analytic solution for rectangular load (Hetenyi 24', 24")
xSol = plateBreakPoint:(1e3):300e3;
y = zeros(size(xSol));
for iX = 1:numel(xSol)
    xNow = abs(xSol(iX)-plateBreakPoint); % distance from breakpoint
    Q = qx(1)/(2*k);
    Alx = exp(-lambda*xNow)*(cos(lambda*xNow) + sin(lambda*xNow));
    Blx = exp(-lambda*xNow)*sin(lambda*xNow);
    Bll = exp(-lambda*l)*sin(lambda*l);
    Cll = exp(-lambda*l)*(cos(lambda*l)-sin(lambda*l));
    Dlx = exp(-lambda*xNow)*cos(lambda*xNow);
   
    if xNow <= l
        % points under load (24')
        Dlm = exp(-lambda*(l-xNow))*cos(lambda*(l-xNow));
        y(iX) = Q*((1+Bll-Cll)*Alx - (1+2*Bll-Cll)*Blx + (2-Dlx-Dlm));

    else 
        % points right of load (24")
        Dlm = exp(-lambda*(xNow-l))*cos(lambda*(xNow-l));
        y(iX) = Q*((1+Bll-Cll)*Alx - (1+2*Bll-Cll)*Blx - (Dlx-Dlm));
    
    
    end
end

% superposition solution
ySuper = flex1d(x,qx,xSol,lambda,k,'broken',plateBreakPoint);

% calculate error
errSuperRectangle = (1/numel(ySuper))*sum((ySuper - y).^2);
fprintf('Mean-squared error: %3.4e (km)\n',errSuperRectangle/1e3);

% calculate airy isostatic deflection:
airyfnc = @(b1) (Tsed-b1)*(rho_sed)/(rho_m-rho_sed)-b1;
airyDef = fzero(airyfnc,0);
% fprintf('Airy deflection: %3.2f (km)\n',airyDef/1e3);

hFig = figure(1); clf
hAx(1) = subplot(3,1,1);
qPlot = [qx,0,0];
xPlot = [x,x(end),xSol(end)];
plot(xPlot/1e3,qPlot/1e6);
ylabel('load (MPa)')
xlabel('x (km)')
title(sprintf('Rectangular load, plate broken at -50 km'))

hAx(2) = subplot(3,1,2);
hold on
plot(xSol/1e3,-ySuper/1e3,'-k','LineWidth',8)
plot(xSol/1e3,-y/1e3,'LineWidth',2)

airyPlotX = [x(1),x(end),x(end),xSol(end)];
airyPlotY = [airyDef,airyDef,0,0];
plot(airyPlotX/1e3,-airyPlotY/1e3,'--k')
ylabel('Deflection (km)')
xlabel('x (km)')
box on
arrow([60,0],[x(end)/1e3,0],'Length',8,'Width',2);
hAiryT = text(-45,0,'Airy Isostacy Depth');
hAiryT.Parent = hAx(2);
hAiryT.HorizontalAlignment = 'left';
publishfigure(hAiryT);

% hAnT = text(hAiryT.Extent(1),-1,'Analytic Solution');
% publishfigure(hAnT);
% hNumT = text(hAiryT.Extent(1),-1.5,'Superposition Solution');
% publishfigure(hNumT);
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
% savefigure_cjt(hFig,'github_benchmark3','-png')