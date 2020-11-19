%clear all

density = 50000;

c = Context(4)
c = cs(c,'p.e',0.0549006)
c = cs(c,'p.f0',0)
c = cs(c,'p.mu',0.012150581623434)
c = solveTrueAnomaly(density,c)
sol = cg(c,'p.fsol')
tspan = linspace(0,2*pi,density);
fhist = deval(sol,tspan);

e = cg(c,'p.e');
mu = cg(c,'p.mu');

for i = 1:size(tspan,2)
    t = tspan(i);
    f = fhist(i);

    A = [cos(t) sin(t)
    -sin(t) cos(t)];

    xm1inertial = - mu * cos(f) / (1 + e*cos(f));
    ym1inertial = - mu * sin(f) / (1 + e*cos(f));
    
    m1posinertial(:,i) = [xm1inertial
                          ym1inertial];
    
    m1pos(:,i) = A*m1posinertial(:,i);

    xm2inertial = (1 - mu) * cos(f) / (1 + e*cos(f));
    ym2inertial = (1 - mu) * sin(f) / (1 + e*cos(f));
    
    m2posinertial(:,i) = [xm2inertial
                          ym2inertial];
    
    m2pos(:,i) = A*m2posinertial(:,i);
end

hold on
scatter(m1posinertial(1,:),m1posinertial(2,:),'.b')
scatter(m2posinertial(1,:),m2posinertial(2,:),'.r')
scatter(m1pos(1,:),m1pos(2,:),'.k')
scatter(m2pos(1,:),m2pos(2,:),'.k')