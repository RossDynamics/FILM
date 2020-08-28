function sphere = threesphere(n)
%Generates a unit 3-sphere (aka a hypersphere) with 
%point density n. I know that this implementation isn't stellar...
%at some point, it probably could be improved
psi = linspace(0, pi, n);
theta = linspace(0, pi, n);
phi = linspace(0, 2*pi, n);

[PSI,THETA,PHI] = ndgrid(psi,theta,phi);

x0 = cos(PSI);
x1 = sin(PSI) .* cos(THETA);
x2 = sin(PSI) .* sin(THETA) .* cos(PHI);
x3 = sin(PSI) .* sin(THETA) .* sin(PHI);

sphere = [x0(:) x1(:) x2(:) x3(:)];
sphere = unique(sphere, 'rows').';

end

