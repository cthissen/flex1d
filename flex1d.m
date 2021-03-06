function [deflection] = flex1d(xLoad,qLoad,xSol,lambda,k,plateCase,xBroken)
%   flex1d(xLoad,qLoad,xSol,lambda,k,plateCase,xBroken) calculates the
%   deflection of a beam on an elastic foundation. xLoad and qLoad must be
%   vectors that define the location and magnitude of the load. xSol is a
%   vector that gives the location where the deflection will be calculated.
%   lambda and k are scalars that define the flexural response. plateCase
%   is a string that must be either 'infinite' or 'broken'. xBroken defines
%   the location of the break for the broken plate solution. deflection
%   gives the flexural deflection of the beam in the same units as xLoad
%   and xSol. 
%
%   [deflection] = flex1d(xLoad,qLoad,xSol,lambda,k,'infinite') uses the infinite beam
%   solutions
%
%   [deflection] = flex1d(xLoad,qLoad,xSol,lambda,k,'broken',xBroken) uses the broken beam
%   solutions, with the break at xBroken.
%
%
%   For example. Calculate the deflection of an infinite beam with a
%   rectangular load:
%     g = 9.81;      % m/s^2, gravity
%     rho_m = 3340;  % kg/m^3, density of material below beam. Determines bouyancy force. 
% 
%     % set flexure parameters
%     lambda = 1/54000; % flexural parameter (1/meter)
%     k = g*rho_m;      % elastic foundation parameter 
% 
%     % define loading
%     Tsed = 3e3;    % m, thickness of sedimentary load
%     rho_sed = 2700;% kg/m^3, density of material that infills deflection (e.g. sediments)
% 
%     % Define loading
%     %... load extends from qxLeft to qxRight with constant magnitude q0
%     dx = 1e2; % meters
%     qxLeft  = -50e3;  % meters
%     qxRight  = 100e3; % meters
%     x = qxLeft:dx:qxRight; % load vector
%     q0 = rho_sed*Tsed*g;
%     qx = q0 + zeros(size(x)); % load, kg/m^2 ? 
% 
%     % calculate deflection
%     xSol = -300e3:(1e3):300e3; % vector of locations where we want to calculate deflection
%     deflection = flex1d(x,qx,xSol,lambda,k,'infinite'); % calculated deflection
% 
%     % plotting
%     hFig = figure(1); clf
%     subplot(2,1,1)
%     qPlot = [0,0,qx,0,0];
%     xPlot = [xSol(1),x(1),x,x(end),xSol(end)];
%     plot(xPlot/1e3,qPlot/1e6);
%     ylabel('load (MPa)')
%     xlabel('x (km)')
%     title(sprintf('Rectangular Load'))
% 
%     subplot(2,1,2)
%     plot(xSol/1e3,-deflection/1e3)
%     ylabel('Deflection (km)')
%     xlabel('x (km)')
% 
%
%   Additional examples can be found in test_example1.m. test_example2.m,
%   test_example3.m, test_benchmark1.m, test_benchmark2.m,
%   test_benchmark3.m
% 
%   Christopher Thissen (cthissen@gmail.com)
%   Yale University
%   Feb 10, 2016

%% check and parse inputs
narginchk(6,7);
validateattributes(xLoad,{'numeric'},{'vector'});
validateattributes(qLoad,{'numeric'},{'vector'});
validateattributes(xSol,{'numeric'},{'vector'});
validateattributes(lambda,{'numeric'},{'scalar'});
validateattributes(k,{'numeric'},{'scalar'});
validatestring(plateCase,{'infinite','broken'});

switch plateCase
    case 'infinite'
        narginchk(6,6);

    case 'broken'
        narginchk(7,7);
        validateattributes(xBroken,{'numeric'},{'scalar'});
end        

checksamesize(xLoad,qLoad);

%% Calculate deflection

switch plateCase
    case 'infinite'

        [deflection] = flexure_infiniteplate(xLoad,qLoad,xSol,lambda,k);
        
        
    case 'broken'


        % first calculate bending moment and shear at xBroken
        [~,MA,QA] = flexure_infiniteplate(xLoad,qLoad,xBroken,lambda,k);
        
        % calculate end-conditioning loads that cancel out load and moment
        % at xBroken
        P0 = 4*(lambda*MA + QA);
        M0 = (-2/lambda)*(2*lambda*MA + QA);
        
        % calculate infinite plate solution
        [dInfinite,~,~] = flexure_infiniteplate(xLoad,qLoad,xSol,lambda,k);
        
        % now calculate addtional deflection due to end conditioning forces
        xNow = abs(xSol - xBroken);
        A = Alamx(lambda,xNow);
        B = Blamx(lambda,xNow);
        brokenEndDefl = (P0*lambda)/(2*k)*A; % deflection due to P0
        brokenEndMomD = (M0*lambda^2/k)*B;   % deflection due to M0

        % total deflection
        deflection = dInfinite + brokenEndDefl + brokenEndMomD; 
end

% validate output
checksamesize(xSol,deflection);



end

function [deflection,moment,shearForce] = flexure_infiniteplate(x,qx,xSol,lambda,k)
% superposition solution (Hetenyi 1949 pg 18)

deflection = zeros(size(xSol));
moment     = deflection;
shearForce = deflection;

for iX = 1:numel(xSol)  % loop over all points
    if numel(x) == 1
        % point solution
        xNow = abs(xSol(iX)-x);
        deflection(iX) = qx*lambda/(2*k) * Alamx(lambda,xNow);
        
    else

        % integrate contributions from all point loads to deflection at xSol(iX)
        xAll = abs(xSol(iX)-x);
        dVal = (lambda/(2*k))*qx.*Alamx(lambda,xAll);
        dVal2 = trapz(x,dVal); % use x here to get the correct spacing
        deflection(iX) =  dVal2;
        
        % also calculate moment and shear
        mVal = (qx/(4*lambda)).*Clamx(lambda,xAll);
        sVal  = (qx/2).*Dlamx(lambda,xAll);

        mVal2 = trapz(x,mVal);
        sVal2 = trapz(x,sVal);
        
        moment(iX) = mVal2;
        shearForce(iX) = sVal2;        
    end
    
end
end

%% Dependent Functions
function [A] = Alamx(lambda,x)
% Hetenyi function A_{\lambda x}
lx = lambda*x;
A = exp(-lx).*(cos(lx)+sin(lx));    
end

function [B] = Blamx(lambda,x)
% Hetenyi function B_{\lambda x}
lx = lambda*x;
B = exp(-lx).*sin(lx);
end

function [C] = Clamx(lambda,x)
% Hetenyi function A_{\lambda x}
lx = lambda*x;
C = exp(-lx).*(cos(lx)-sin(lx));
    
end

function [D] = Dlamx(lambda,x)
% Hetenyi function D_{\lambda x}
lx = lambda*x;
D = exp(-lx).*cos(lx);
    
end

function [] = checksamesize(data1,data2)
% check if data1 and data2 have the same dimensions

%%
narginchk(2,2);
validateattributes(data1,{'numeric'},{});
validateattributes(data2,{'numeric'},{});

%%
nData1 = size(data1);
nData2 = size(data2);

nDims1 = numel(nData1);
nDims2 = numel(nData2);

if nDims1 ~= nDims2
    error('data do not have the same number of dimensions');
end
if sum(nData1 == nData2) ~= nDims1
    error('data have the same number of dimensions but are not the same size');
end

end