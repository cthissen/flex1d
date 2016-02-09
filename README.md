# flex1d

What is it?
----------------- 
flex1d calculates the flexure of the lithosphere and other elastic plates under arbitrary loading.

The code uses solutions from Hetényi, M. (1971). Beams on elastic foundation: theory with applications in the fields of civil and mechanical engineering. University of Michigan.

Applications to the lithosphere are described in Watts, A. B. (2001). Isostasy and Flexure of the Lithosphere. Cambridge University Press.

For comments, questions, or suggestions, please email cthissen@gmail.com or 
leave a comment under the issues tab at github.com/cthissen/estream2

Christopher J. Thissen, Yale University  

Usage 
-----------------
The syntax is simple:
````
% Rectangular distributed load on an infinite plate

% Define parameters
g = 9.81;       % gravity (meters/seconds^2)
rho_m = 10;     % density of material below beam. Determines bouyancy force.
rho_infill = 1; % density of material that infills deflection (e.g. water or sediments)
lambda = 1/54;  % flexural parameter. (1/meters)
k = g*(rho_m - rho_infill); % "elastic foundation" parameter. For the lithosphere this is the bouyancy force.

% Define loading
dx = 1; % meters
qxLeft  = -50; % meters
qxRight  = 150; % meters
x = qxLeft:dx:qxRight; % load vector
qx = 1e3 + zeros(size(x)); % load parameters
xSol = -300:(1*dx):300; % vector where we want the flexure

% calculate deflection
deflection = flex1d(x,qx,xSol,lambda,k,'infinite');

figure(1); clf
plot(xSol,-deflection)
````

Additional examples can be found in test_flex1d.m

Calculating Lithospheric Flexure Parameters
-----------
The code takes as input the somewhat cryptic parameters lambda and k. 

### Calculating lambda
Lambda describes the flexural response of the lithosphere. A typical value is 1/54000, in units of meters. 

We usually have estimates for the elastic thickness, Te. This can be converted the lambda by first converting to the flexure rigidity, 

<img>

where E is Young's modulous (typically 10e10 Pascals) and v is Poisson's ratio (0.25). 

Lambda is then related to D by


### Calculating k
The solutions are for a beam on an elastic foundation. This means that downward deflection of the plate is resisted both by the flexural strength of the beam itself, and the force to move the underlying material out of the way. The elastic part of the foundation means that, like a spring, the resisting force from the foundation is a linear function of the distance, F = kx. For lithospheric problems, the resisting force of the foundation is due to buoyancy. As the deflection pushes the mantle out of the way, the difference in densities creates a bouyance force that resists the deflection. So how to calculate k?

<eqn>

The relevant densities are the density of the mantle, and the density of the material that infills the deflection, such as water or sediments. 

<img>

Make sure your units match! If you use g = 9.81 (m/s^2), then your densities should be kg/m^3. 

### Example
Here's an example that calculates lambda and k using typical lithospheric values.
````


````


# Infinite vs Broken Plates
The code will provide solutions for both infinite and broken plates. Infinite plates extend towards infinity in both directions. Broken plates are "broken" at a specific horiztonal position. Broken here means that the plate has zero moment and shear where it is broken. Other types of end conditions (e.g. hinged, or fixed location) are also possible. If there's interest these solutions can be easily added. Just send me an email (cthissen@gmail.com) or leave a comment on the issues tab above. 

You must specify whether your plate is infinite or broken.  
Here's the syntax for an infinite plate
````
deflection = flex1d(x,qx,xSol,lambda,k,'infinite');
````
And here's the syntax for a broken plate
````
deflection = flex1d(x,qx,xSol,lambda,k,'broken',plateBreakPoint);
````
where plateBreakPoint is the x location (in the same units) of the break in the plate. 

Here's the first example, but using a plate broken at the left-limit of the load.
````
% Rectangular distributed load on a broken plate

% Define parameters
g = 9.81;       % gravity (meters/seconds^2)
rho_m = 10;     % density of material below beam. Determines bouyancy force.
rho_infill = 1; % density of material that infills deflection (e.g. water or sediments)
lambda = 1/54;  % flexural parameter. (1/meters)
k = g*(rho_m - rho_infill); % "elastic foundation" parameter. For the lithosphere this is the bouyancy force.

% Define loading
dx = 1; % meters
qxLeft  = -50; % meters
qxRight  = 150; % meters
x = qxLeft:dx:qxRight; % load vector
qx = 1e3 + zeros(size(x)); % load parameters
xSol = qxLeft:(1*dx):300; % vector where we want the flexure

% calculate deflection
deflection = flex1d(x,qx,xSol,lambda,k,'broken',qxLeft);

figure(1); clf
plot(xSol,-deflection)
````
