%Used for analyzing trajectories with components on the center manifold
%using the symplectic eigenbasis.
%Run analyzeMonodromy first.

%Let's set the eigenbasis first
c = coordset(c,eigenbasis);

q1coord = 1e-7;
p1coord = 1e-5;

circSize = 1e-4;

n = 500;

circ = circSize * onesphere(n);

coords = [q1coord * ones(n,1) circ(:,1) p1coord * ones(n,1) circ(:,2)].'

c = cs(c,'s.o.v.dmode','2');
close all;
cplot(coords,c,'ob');
c = cs(c,'s.o.v.dmode','1');

coordget(c)

c = cs(c,'lm.tube',coords);

c = coordreset(c);

c = useMomentum(c);

coordget(c)

tubecoords = cg(c,'lm.tube')

%Now, we integrate. For speed, we enable caching.
c = startCaching(c);

periods = 1;
numPts = 2;

tspan = linspace(0,periods * getSolarPeriod(c),numPts);

figure
hold on 

fdcoords = [];
bdcoords = [];

for i = 1:size(tubecoords,2)
    [sol,c] = integ(tspan,tubecoords(:,i),c);
    y = deval(sol,tspan);
    p = cplot(y,c,'bo');
    p.Color(4) = 0.25;
    
    fdcoords = [fdcoords y(:,end)];
end

%fdshp1.Alpha = Inf;
%fdshp2 = alphaShape(fdcoords(2,:)',fdcoords(4,:)');
%fdshp2.Alpha = Inf;
%fdArea = 
%plot(fdshp1);

for i = 1:size(tubecoords,2)
    [sol,c] = integ(-tspan,tubecoords(:,i),c);
    y = deval(sol,-tspan);
    p = cplot(y,c,'bo');
    p.Color(4) = 0.25;
    
    bdcoords = [bdcoords y(:,end)];
end

%bdshp1 = alphaShape(bdcoords(1,:)',bdcoords(3,:)');
%bdshp1.Alpha = Inf;
%plot(bdshp1);

icInvariant = - polyarea(tubecoords(1,:)',tubecoords(3,:)') + ...
               polyarea(tubecoords(2,:)',tubecoords(4,:)')
fdInvariant = - polyarea(fdcoords(1,:)',fdcoords(3,:)') + ...
               polyarea(fdcoords(2,:)',fdcoords(4,:)')
bdInvariant = polyarea(bdcoords(1,:)',bdcoords(3,:)') - ...
               polyarea(bdcoords(2,:)',bdcoords(4,:)')

c = stopCaching(c);