function U = BCM_U(x, y, t, mu, mu0, a0, omegam0, thetam00)
%BCM_U The potential energy in the BCM expressed in velocity coordinates.

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
r1 = sqrt((x+mu2)^2 + y^2);
r2 = sqrt((x-mu1)^2 + y^2);

%From where does this extra potential term come? Is it because the body is
%moving? See the Hamiltonian in Simó et al. (1995).
U = mu0 * (x * cos(thetam0) + y * sin(thetam0)) / a0^2 ...
                           - mu0 / r0 - mu1 / r1 - mu2 / r2; 

end

