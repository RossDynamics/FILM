function U = ER3BP_U(x, y, t, mu, f, e)
%ER3BP_U The potential energy in the ER3BP expressed in velocity 
%coordinates. f should be a function.

%--Variables--

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

%Should this potential have extra terms as in the BCM?
U = - (1 - mu) / r1 - mu / r2; 

end

