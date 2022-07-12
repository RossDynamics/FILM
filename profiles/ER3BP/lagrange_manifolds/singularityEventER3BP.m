function [value,isterminal,direction] = singularityEventER3BP(t,y,mu,ac,e)
%SINGULARITYEVENT An event function that prevents hasflm from integrating
%too close to a mass-based singularity. An active coordinates struct must
%be passed to singularity event.

%--Variables--

f = y(5);

%The rotation matrix from the inertial frame to the rotating frame
A = [cos(t) sin(t)
    -sin(t) cos(t)];

%m1 positioning
xm1inertial = - mu * cos(f) / (1 + e*cos(f));
ym1inertial = - mu * sin(f) / (1 + e*cos(f));

m1 = A*[xm1inertial
           ym1inertial];
       
%m2 positioning
xm2inertial = (1 - mu) * cos(f) / (1 + e*cos(f));
ym2inertial = (1 - mu) * sin(f) / (1 + e*cos(f));

m2 = A*[xm2inertial
           ym2inertial];

%We have to get the location of y in standard coordinates
ysta = new2standard(y,ac.basis.value,ac.origin.value);

%We create circles around each singularity that the integrator should not
%cross.
value(1) = norm(ysta(1:2) - m1) - 5e-2;
isterminal(1) = true;
direction(1) = 0;

value(2) = norm(ysta(1:2) - m2) - 5e-2;
isterminal(2) = true;
direction(2) = 0;

end

