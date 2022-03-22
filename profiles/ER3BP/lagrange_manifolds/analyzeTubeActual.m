%Used for analyzing trajectories with components on the center manifold
%using the symplectic eigenbasis.
%Run analyzeMonodromyER3BP and analyzeTransitER3BP first.

%Let's set the eigenbasis first
c = coordset(c,eigenbasis);

%We get the transit initial conditions and then create a circle with radius
%equal to the distance of the default center manifold points to the origin
transit = cg(c,'lm.transit');

coords = [];
for i = 1:size(transit,2)
    q1coord = transit(1,i);
    p1coord = transit(3,i);
    
    circSize = sqrt(transit(2,i)^2+transit(4,i)^2);
    
    % q1coord = 1e-7;
    % p1coord = 1e-5;
    %
    % circSize = 1e-4;
    
    n = 80;
    
    circ = circSize * onesphere(n);
    
    coords = [coords [q1coord * ones(n,1) circ(:,1) p1coord * ones(n,1) circ(:,2) zeros(n,1)].'];
end

c = cs(c,'s.o.v.dmode','2');
close all;
cplot(coords,c,'.b');
c = cs(c,'s.o.v.dmode','position');

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

tspan = linspace(0,periods * 2*pi,numPts);

figure
hold on 

fdcoords = [];
bdcoords = [];

for i = 1:size(tubecoords,2)
    [solfd,c] = integ(tspan,tubecoords(:,i),c);
    yfd = deval(solfd,tspan);
    fdcoords = [fdcoords yfd(:,end)];
end

p = cplot(fdcoords,c,'b.');
p.Color(4) = 0.25;

set(gca,'FontSize',20)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
axis square
axis equal

figure;

%fdshp1.Alpha = Inf;
%fdshp2 = alphaShape(fdcoords(2,:)',fdcoords(4,:)');
%fdshp2.Alpha = Inf;
%fdArea = 
%plot(fdshp1);

for i = 1:size(tubecoords,2)
    [solbd,c] = integ(-tspan,tubecoords(:,i),c);
    ybd = deval(solbd,-tspan);
    bdcoords = [bdcoords ybd(:,end)];
end

p = cplot(bdcoords,c,'b.');
p.Color(4) = 0.25;

set(gca,'FontSize',20)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
axis square
axis equal

%bdshp1 = alphaShape(bdcoords(1,:)',bdcoords(3,:)');
%bdshp1.Alpha = Inf;
%plot(bdshp1);

% icInvariant = - polyarea(tubecoords(1,:)',tubecoords(3,:)') + ...
%                polyarea(tubecoords(2,:)',tubecoords(4,:)')
% fdInvariant = - polyarea(fdcoords(1,:)',fdcoords(3,:)') + ...
%                polyarea(fdcoords(2,:)',fdcoords(4,:)')
% bdInvariant = polyarea(bdcoords(1,:)',bdcoords(3,:)') - ...
%                polyarea(bdcoords(2,:)',bdcoords(4,:)')

c = stopCaching(c);