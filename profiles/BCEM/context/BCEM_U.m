function U = BCEM_U(x, y, t, mu, f, e, mu0, a0, omegam0, thetam00)
%ER3BP_U The potential energy in the BCEM expressed in velocity 
%coordinates. f should be a function.

%---BCM Contribution---

%Masses
mu1 = 1 - mu;
mu2 = mu;
    
%--Variables--

%Solar positioning
thetam0 = -omegam0 * t + thetam00;
    
xm0 = a0 * cos(thetam0);
ym0 = a0 * sin(thetam0);

%Particle state
r0 = sqrt((x-xm0)^2 + (y-ym0)^2);

%---ER3BP Contribution---

%The rotation matrix from the inertial frame to the rotating frame
A = [cos(t) sin(t)
    -sin(t) cos(t)];

%m1 positioning
xm1inertial = - mu * cos(f(t)) / (1 + e*cos(f(t)));
ym1inertial = - mu * sin(f(t)) / (1 + e*cos(f(t)));

m1pos = A*[xm1inertial
           ym1inertial];
       
%m2 positioning
xm2inertial = (1 - mu) * cos(f(t)) / (1 + e*cos(f(t)));
ym2inertial = (1 - mu) * sin(f(t)) / (1 + e*cos(f(t)));

m2pos = A*[xm2inertial
           ym2inertial];

%Particle state
r1 = sqrt((x-m1pos(1))^2 + (y-m1pos(2))^2);
r2 = sqrt((x-m2pos(1))^2 + (y-m2pos(2))^2);

%---Potential---
U = mu0 * (x * cos(thetam0) + y * sin(thetam0)) / a0^2 ...
                           - mu0 / r0 - mu1 / r1 - mu2 / r2;

end

