function U = CR3BP_U(x, y, t, mu)
%BCM_U The potential energy in the CR3BP expressed in velocity coordinates.

%Masses
mu1 = 1 - mu;
mu2 = mu;
    
%--Variables--
    
%Particle state
r1 = sqrt((x+mu2)^2 + y^2);
r2 = sqrt((x-mu1)^2 + y^2);

%From where does this extra potential term come? Is it because the body is
%moving? See the Hamiltonian in Simó et al. (1995).
U = - mu1 / r1 - mu2 / r2; 

end

