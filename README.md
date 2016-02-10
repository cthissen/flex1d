# flex1d

What is it?
----------------- 
flex1d calculates the flexure of the lithosphere and other elastic plates under arbitrary loading. The solution is implemented in 1 dimension, so it calculates the deflection for a beam of any length and unit width. 

The code uses solutions from Hetényi, M. (1971). Beams on elastic foundation: theory with applications in the fields of civil and mechanical engineering. University of Michigan.

Applications to the lithosphere are described in Watts, A. B. (2001). Isostasy and Flexure of the Lithosphere. Cambridge University Press.

For comments, questions, or suggestions, please email cthissen@gmail.com or 
leave a comment under the issues tab at github.com/cthissen/estream2

Christopher J. Thissen, Yale University  

Usage 
-----------------
The syntax is simple. For an infinite plate:  
````
deflection = flex1d(x,qx,xSol,lambda,k,'infinite');
````  
For a broken plate:  
````
deflection = flex1d(x,qx,xSol,lambda,k,'broken',plateBreakPoint);
````
x gives the location of the load   
qx gives the value of the load  
xSol is the location where we want to calculate the deflection  
lambda describes the flexural response of the plate  
k describe the resistance of the underlying material to deflection  
'infinite' or 'broken' specifies whether the plate extends to infinity or is broken.  
plateBreakPoint is the x location of the break in the plate.   
deflection is the deflection of the beam

#### Infinite vs Broken Plates
The code will provide solutions for both infinite and broken plates. Infinite plates extend towards infinity in both directions. Broken plates are "broken" at a specific horiztonal position. Broken means that the plate has zero moment and shear where it is broken. Other types of end conditions (e.g. hinged, or fixed location) are also possible. If there's interest these solutions can be easily added. Just send me an email (cthissen@gmail.com) or leave a comment on the issues tab above. 

Calculating Lithospheric Flexure Parameters
-----------
The code takes as input the somewhat cryptic parameters lambda and k. 

Lambda describes the flexural response of the lithosphere. A typical value is 1/54000, in units of meters. 

For the lithosphere we usually have estimates for the elastic thickness, Te, not lambda. Te can be converted the lambda by first converting to the flexure rigidity,   
<img src="https://github.com/cthissen/flex1d/blob/master/images/D.png" alt="alt text" height="50px">  
where E is Young's modulous (typically 100 Gigapascals) and v is Poisson's ratio (0.25). 

Lambda is related to D by  
<img src="https://github.com/cthissen/flex1d/blob/master/images/lambda.png" alt="alt text" height="50px"> 

The downward deflection of the plate is resisted both by the flexural strength of the beam itself, and the force needed to move the underlying material out of the way. For lithospheric problems, the resisting force of the foundation is due to buoyancy. As the deflection pushes the mantle out of the way, the difference in densities creates a bouyance force that resists the deflection. K is calculated simply as  

<img src="https://github.com/cthissen/flex1d/blob/master/images/k.png" alt="alt text" height="25px"> 

where <img src="https://github.com/cthissen/flex1d/blob/master/images/rhom.png" alt="alt text" height="20px"> is the density of the mantle. Make sure your units match! If you use g = 9.81 (m/s^2), then your densities should be kg/m^3. 



### Example
Here's an example that calculates lambda and k using typical lithospheric values.  
<img src="https://github.com/cthissen/flex1d/blob/master/images/github_ex2.png" alt="alt text" width="400px"> 

The above image is created using this code:  
````
%% Infinite beam, distributed load
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
xSol = -300e3:(1e3):300e3; % vector of locations where we want to calculate deflection
deflection = flex1d(x,qx,xSol,lambda,k,'infinite'); % calculated deflection

% plotting
hFig = figure(1); clf
subplot(2,1,1)
qPlot = [0,0,qx,0,0];
xPlot = [xSol(1),x(1),x,x(end),xSol(end)];
plot(xPlot/1e3,qPlot/1e6);
ylabel('load (MPa)')
xlabel('x (km)')
title(sprintf('Rectangular Load'))

subplot(2,1,2)
plot(xSol/1e3,-deflection/1e3)
ylabel('Deflection (km)')
xlabel('x (km)')
````

Here's the same example, but using a plate broken at -50 km.  
<img src="https://github.com/cthissen/flex1d/blob/master/images/github_ex3.png" alt="alt text" width="400px"> 

````
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
title(sprintf('Rectangular Load'))

subplot(2,1,2)
plot(xSol/1e3,-deflection/1e3)
ylabel('Deflection (km)')
xlabel('x (km)')
````
#### A note about the density of the infilling material


## Benchmarks
Hetenyi (1949) provides a number of analytical solutions that can be used to test the implementation. 


