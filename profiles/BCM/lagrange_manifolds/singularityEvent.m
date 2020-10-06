function [value,isterminal,direction] = singularityEvent(t,y,mu,ac)
%SINGULARITYEVENT An event function that prevents hasflm from integrating
%too close to a mass-based singularity. An active coordinates struct must
%be passed to singularity event.

m1 = [-mu 0].';
m2 = [1 - mu 0].';

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

