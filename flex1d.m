function [deflection] = flex1d(xLoad,qLoad,xSol,lambda,k,plateCase,xBroken)
% calculate the felxural response of an elastic plate on an elastic
% foundation. 
% Assumes
% 1. Plate is perfectly elastic, with no viscous deformation or fracturing.
% 2. Linear elastic rheology.
% 3. Constant elastic properties
% 4. Underlying material reponds with a force that is linearly related to
% the deflection (Winkler foundation). In the mantle, this is provided by
% the bouyancy force. 
% 5. Isotropic plate, which allows the elasticity to captured by two
% parameters

% Plasti uses k = g*rho_m, also leaves rhof out of alpha
%
% c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% C  NOTE: there is no rhof in the def of alpha.
% c 	see project notes for description of method of calculating flexure and
% c	why the rhof is left out. In brief, it is left out because the forces
% c	acting on the imaginary plate are calculated as the load from the crust
% c	and the load from the water.  Another way to do this problem would be 
% c	to use rhof in the eqn and calculate just the loads from the crust.
% c	In this case the force from any portion of a colm of crust that is below a 
% c	defined sea level is (rhoc-rhof)*g*h and the force form the portion of 
% c	the same colm above sea level (if there is a sub aerial portion) is 
% c	rhoc*g*h', where h' is the height of the colm above sea level
% c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% c calculate flexural parameters
%       alpha1=(4.0*prig/((rhom)*g))**0.25
%       alpha2=(4.0*rrig/((rhom)*g))**0.25
%       plam1=1.0/alpha1
%       plam2=1.0/alpha2
%       fk=rhom*g

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
        % point solution, assumes line load (not qx*dx)
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