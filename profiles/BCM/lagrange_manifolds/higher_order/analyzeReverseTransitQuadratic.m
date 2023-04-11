%Used for analyzing transit using the symplectic eigenbasis, but going
%backwards. Run analyzeMonodromy or analyzeMonodromyER3BP first--
%this code is designed to work with both the BCM and the ER3BP.

%model = "BCM"
model = "ER3BP"

%Change depending on the model:
switch model
    case "BCM"
        period = getSolarPeriod(c) 
        conststart = mp(1e-7);
        constend = mp(4.9e-6);
        h = mp(1e-13);
        k = mp(1e-10);
        
    case "ER3BP"
        period = 2*pi
        conststart = mp(1e-8);
        constend = mp(5e-7);
        h = mp(1e-14);
        k = mp(1e-7);
end

%To consider initial conditions in the upper halfplane, set halfplane equal
%to 1. To consider initial conditions in the lower halfplane, set
%halfplane equal to -1.
halfplane = -1;

%We create an array of initial conditions from within the
%symplectic eigenbasis along the bounding line. We set const (where p1 =
%const + q1) small and h = H2(x) small.

%Explanation of k from the old use of the variable:
%We also set a beginning for the array's q1's. In the linear case, we 
%typically take this beginning to coincide with the p1 + q1 = 0 line, 
%but we want the ability to zoom in, so we take a parameter k such that we
%will start at p1 + k * q1 = 0 and we will end at p1 - k * q1 = 0. As a
%result, set k = 1 for the standard case.

%Current explanation of k:
%All p1's are generated from -k to k.

c = coordreset(c);
c = useMomentum(c);

%We get the quadratic map. It is faster if we try to compute it within the
% %standard momentum frame and then include conversion functionality manually
% disp('Computing phi')
% c = cs(c,'s.i.odeopts',odeset('OutputFcn',@odeprint,'RelTol',...
%                               3e-14,'AbsTol',1e-15));
% [phi,~,~] = stm(linspace(0,-period,10),mp(cg(c,'lm.y0')),c,[],2);
% c = cs(c,'s.i.odeopts',odeset('OutputFcn',[],'RelTol',...
%                               3e-14,'AbsTol',1e-15));
% disp('Done computing phi')

eigp = eigenbasis.basis.value;
eigo = eigenbasis.origin.value;

quadratstd = @(x)(phi{1}(1:4,1:4,end)*x +...
    getfield(1/2*ttv(ttv(tensor(phi{2}(1:4,1:4,1:4,end)),x,2),x,2),'data'));
quadrat=@(x) eigp \ quadratstd(eigp * x);

linstd = @(x)phi{1}(1:4,1:4,end)*x;
lin=@(x) eigp \ linstd(eigp * x);

%Let's set the eigenbasis now
c = coordset(c,eigenbasis);

%Old:
%The other end of the array will not be found dynamically as in the
%standard analyzer. Instead, it will occur when the p1 - k * q1 = 0 line is
%reached; in other words, q1 = const/(k-1). We don't need to use an array
%step but rather we can just preallocate the q1's based on a fixed number 
%of points to be analyzed.

numConsts = 30;
numPoints = 100;
% numConsts = 10;
% numPoints = 10;

constarray = linspace(halfplane*conststart,halfplane*constend,numConsts)';
%constarray = halfplane.*logspace(-10,-7,numConsts)';

arrayics = NaN(getnExtended(c),numPoints*numConsts,'mp');
arrayicsfd = NaN(getnExtended(c),numPoints*numConsts,'mp');
arrayicsfdlin = NaN(getnExtended(c),numPoints*numConsts,'mp');

for j = 1:numConsts
    for i = 1:numPoints
        
        
        %q1array = linspace(-constarray(j) / (k + 1),constarray(j)/(k-1),numPoints);
        p1array = linspace(-k,k,numPoints);
        
        ic = ics_energy_boundary(p1array(i),constarray(j),h,...
                                   cg(c,'p.sigma'),cg(c,'p.a'),...
                                   cg(c,'p.T'),true);
        %If we get complex values, we know we've reached the boundary,
        %and so we break and leave all remaining entries (including this one)
        %as NaN's
        if ~isreal(ic)
            %break;
        end
        
        arrayics(:,i + (j - 1)*numPoints) = [ic
                                             zeros(getnExtended(c)-4,1)];
        arrayicsfd(:,i + (j - 1)*numPoints) = [quadrat(ic)
                                               zeros(getnExtended(c)-4,1)];
        arrayicsfdlin(:,i + (j - 1)*numPoints) = [lin(ic)
                                               zeros(getnExtended(c)-4,1)];
    end
end

disp('------------')
disp('CONDITION 2:')
disp('Last element of the initial conditions array''s q1*p1')
disp(arrayics(1,end)*arrayics(3,end));

disp('h/lambdatilde:')
disp(h / (1/period*log(cg(c,'p.sigma'))));
disp('------------')

%We use eps instead of 0 because numbers that are almost zero may not be
%treated as 0. We have to multiply by halfplane since the sign of transit
%and nontransit initial conditions depends upon the current halfplane.
%We determine whether an orbit is transit or nontransit based on the side
%of the halfplane to which it deflects under one iterate of the quadratic
%map; one iterate seems to be sufficient for the BCM as the phase space
%divergence rates are very high. We also have to determine whether to use
%the q_1 > 0 / q_1 < 0 or p_1 > 0 / p_1 < 0 halfplanes, based on whether
%we're going forwards or backwards.
hp = 3;
c = cs(c,'lm.nontransit',arrayics(:,arrayicsfd(hp,:) * halfplane < -eps('mp')));
c = cs(c,'lm.transit',arrayics(:,arrayicsfd(hp,:) * halfplane > eps('mp')));
c = cs(c,'lm.manifold',arrayics(:,abs(arrayicsfd(hp,:)) < eps('mp')));

c = cs(c,'lm.nontransitfd',arrayicsfd(:,arrayicsfd(hp,:) * halfplane < -eps('mp')));
c = cs(c,'lm.transitfd',arrayicsfd(:,arrayicsfd(hp,:) * halfplane > eps('mp')));
c = cs(c,'lm.manifoldfd',arrayicsfd(:,abs(arrayicsfd(hp,:)) < eps('mp')));

c = cs(c,'lm.nontransitfdlin',arrayicsfdlin(:,arrayicsfd(hp,:) * halfplane < -eps('mp')));
c = cs(c,'lm.transitfdlin',arrayicsfdlin(:,arrayicsfd(hp,:) * halfplane > eps('mp')));
c = cs(c,'lm.manifoldfdlin',arrayicsfdlin(:,abs(arrayicsfd(hp,:)) < eps('mp')));

%Let's plot the initial conditions in the symplectic eigenbasis's saddle
%plane

close all;

%disp('Running analyzeParabola...')
%analyzeParabola
%disp('analyzeParabola complete.')

lambdatil = log(cg(c,'p.sigma'))/period;
fplot(@(q_1)h/(lambdatil*q_1))

c = cs(c,'s.o.v.dmode','1');

hold on
limScale = 0.1;
line([0 0], [-limScale limScale]);
line([-limScale limScale], [0 0]);
cplot(double(cg(c,'lm.nontransit')),c,'.b');
cplot(double(cg(c,'lm.transit')),c,'.r');
cplot(arrayicsfd(:,1:size(cg(c,'lm.nontransit'),2)),c,'ob');
cplot(arrayicsfd(:,(size(cg(c,'lm.nontransit'),2)+1):end),c,'or');
%cplot(cg(c,'lm.manifold'),c,'ob');

set(gca,'FontSize',20)
xlabel('$q_1$','interpreter','latex')
ylabel('$p_1$','interpreter','latex')
axis([-0.1 0.75 -0.1 0.75]*1e-3)

numPeriods = 1;
timescale = -numPeriods * period;
numPts = 5;

% nontransitics = cg(c,'lm.nontransit');
% transitics = cg(c,'lm.transit');
% manifold = cg(c,'lm.manifold');

c = cs(c,'s.o.v.dmode','position');

%We have to reset the coordinate system for now.
c = coordreset(c);
c = useMomentum(c);

%Now, we integrate. For speed, we enable caching.
c = startCaching(c);

timescale = -2 * period;
numPts = 1000;

nontransitics = double(cg(c,'lm.nontransit'));
transitics = double(cg(c,'lm.transit'));
manifold = double(cg(c,'lm.manifold'));

%pause

integfig = figure
hold on
for i = 1:size(nontransitics,2)
    [p,c] = integplot(linspace(0,timescale,numPts),...
                      nontransitics(:,i),c,'b');
    p.Color(4) = 0.75;
    [p,c] = integplot(linspace(0,-timescale,numPts),...
                      nontransitics(:,i),c,'b');
    p.Color(4) = 0.75;
    disp(i)
    drawnow;
end
for i = 1:size(transitics,2)
    [p,c] = integplot(linspace(0,timescale,numPts),...
                      transitics(:,i),c,'r');
    p.Color(4) = 0.75;
    [p,c] = integplot(linspace(0,-timescale,numPts),...
                      transitics(:,i),c,'r');
    p.Color(4) = 0.75;
    disp(i)
    drawnow
end
cplot(cg(c,'lm.nontransitfd'),c,'ob');
cplot(cg(c,'lm.transitfd'),c,'or');
cplot(cg(c,'lm.nontransitfdlin'),c,'+b');
cplot(cg(c,'lm.transitfdlin'),c,'+r');

% for i = 1:size(manifold,2)
%     [~,c] = integplot(linspace(0,timescale,numPts),...
%                       manifold(:,i),c,'m');
%     [~,c] = integplot(linspace(0,-timescale,numPts),...
%                       manifold(:,i),c,'m');
% end

set(gca,'FontSize',20)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
axis equal

disp('CONDITION 1:')
disp('Assess to see whether it holds visually.')
disp('------------')

c = stopCaching(c);
