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
The only code needed is [flex1d.m](https://github.com/cthissen/flex1d/blob/master/flex1d.m). The files that start with "test_" provide specific examples. 

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



### Examples 
Here's an example that calculates lambda and k using typical lithospheric values.  [Code for this example](https://github.com/cthissen/flex1d/blob/master/test_example2.m)  
<img src="https://github.com/cthissen/flex1d/blob/master/images/github_ex2.png" alt="alt text" width="400px"> 

Here's the same example, but using a plate broken at -50 km.  [Code for this example](https://github.com/cthissen/flex1d/blob/master/test_example3.m)  


<img src="https://github.com/cthissen/flex1d/blob/master/images/github_ex3.png" alt="alt text" width="400px"> 

### [Additional Examples](https://github.com/cthissen/flex1d/wiki/Benchmarks)



