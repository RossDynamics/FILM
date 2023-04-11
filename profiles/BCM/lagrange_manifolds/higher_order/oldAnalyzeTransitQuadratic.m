%Used for analyzing transit using the symplectic eigenbasis.
%Run analyzeMonodromy first. Doesn't work very well and has been superceded
%by the simpler approach in analyzeTransitQuadratic.

%Let's set the eigenbasis first
c = coordset(c,eigenbasis);

%To consider initial conditions in the upper halfplane, set halfplane equal
%to 1. To consider initial conditions in the lower halfplane, set
%halfplane equal to -1.
halfplane = 1;

%We create an array of initial conditions from within the
%symplectic eigenbasis along the bounding line. We set const (where p1 =
%const + q1) small and h = H2(x) small.

const = halfplane * 1e-6;
h = 1e-13;
%h = 0;

%We also set a beginning for the array's q1's. In the linear case, we 
%typically take this beginning to coincide with the p1 + q1 = 0 line, 
%but we want the ability to zoom in, so we take a parameter k such that we
%will start at p1 + k * q1 = 0 and we will end at p1 - k * q1 = 0. As a
%result, set k = 1 for the standard case.

%k = 8;
k = 5e+5;

%The other end of the array will not be found dynamically as in the
%standard analyzer. Instead, it will occur when the p1 - k * q1 = 0 line is
%reached; in other words, q1 = const/(k-1). We don't need to use an array
%step but rather we can just preallocate the q1's based on a fixed number 
%of points to be analyzed.

numPoints = 100;
q1array = linspace(-const / (k + 1),const/(k-1),numPoints);

arrayics = NaN(getnExtended(c),numPoints);

for i = 1:numPoints
    
    ic = ics_energy_boundary(q1array(i),const,h,...
                               cg(c,'p.sigma'),cg(c,'p.a'),cg(c,'p.T'));
    %If we get complex values, we know we've reached the boundary,
    %and so we break and leave all remaining entries (including this one)
    %as NaN's
    if ~isreal(ic)
        break;
    end
    arrayics(:,i) = ic;
end

disp('------------')
disp('CONDITION 2:')
disp('Last element of the initial conditions array''s q1*p1')
disp(arrayics(1,end)*arrayics(3,end));

disp('h/lambdatilde:')
disp(h / (1/getSolarPeriod(c)*log(cg(c,'p.sigma'))));
disp('------------')

%We use eps instead of 0 because numbers that are almost zero may not be
%treated as 0. We have to multiply by halfplane since the sign of transit
%and nontransit initial conditions depends upon the current halfplane.
c = cs(c,'lm.nontransit',arrayics(:,arrayics(1,:) * halfplane < -eps));
c = cs(c,'lm.transit',arrayics(:,arrayics(1,:) * halfplane > eps));
c = cs(c,'lm.manifold',arrayics(:,abs(arrayics(1,:)) < eps));

%Let's plot the initial conditions in the symplectic eigenbasis's saddle
%plane

close all;

%disp('Running analyzeParabola...')
%analyzeParabola
%disp('analyzeParabola complete.')

lambdatil = log(cg(c,'p.sigma'))/getSolarPeriod(c);
fplot(@(q_1)h/(lambdatil*q_1))

c = cs(c,'s.o.v.dmode','1');

hold on
limScale = 0.1;
line([0 0], [-limScale limScale]);
line([-limScale limScale], [0 0]);
cplot(cg(c,'lm.nontransit'),c,'.b');
cplot(cg(c,'lm.transit'),c,'.r');
%cplot(cg(c,'lm.manifold'),c,'ob');

set(gca,'FontSize',20)
xlabel('$q_1$','interpreter','latex')
ylabel('$p_1$','interpreter','latex')
axis([-0.1 0.75 -0.1 0.75]*1e-3)

periods = 1;
timescale = periods * getSolarPeriod(c);
numPts = 5;

nontransitics = cg(c,'lm.nontransit');
transitics = cg(c,'lm.transit');
manifold = cg(c,'lm.manifold');

%Now, we plot the quadratic stable manifold.
%We integrate backwards to get the map yielding the stable manifold
[revphi,~,~] = stm(linspace(0,-getSolarPeriod(c),3),cg(c,'lm.y0'),c,[],2)

revquadrat = @(x)(revphi{1}(1:4,1:4,end)*x +...
    double(1/2*ttv(ttv(tensor(revphi{2}(1:4,1:4,1:4,end)),x,2),x,2)));

%Approximates the quadratic stable manifold. I know I shouldn't be copying
%and pasting code when I could be writing it more modularly, but I'm at the
%end of my dissertation studies, ok?
pointDensity = 50;
epsilon = 1e-14;
%iterates confusingly includes the first iterate, which isn't iterated
iterates = 2;

xquad = [];

nutilde = acos(cg(c,'p.a')) / getSolarPeriod(c);
%centerradguess = sqrt(2 * const / nutilde);
% q2guess = 8.5930e-09
% p2guess = 1.7588e-09
q2guess = -3.031e-9
p2guess = -5.73e-10
% q2guess = 5.55e-10;
% p2guess = 1.047e-10;
% q2guess = 0.8572e-8;
% p2guess = 1.77e-9;
%investigate whether there are multiple sets of initial conditions at same
%energy but different magnitudes of the end phase in the center projection


%p1vals = epsilon * linspace(0,1,pointDensity);
p1vals = epsilon * linspace(0.2,0.23,pointDensity);

%We prepare to refine the center projection coordinates necessary for
%targeting the correct energy
ptMaker = @(p1val,q2,p2)[0
                         q2
                         p1val
                         p2
                         zeros(getnExtended(c)-4,1)];

refineHandle = @(p1val,q2,p2)(H2energy(0,...
                 revquadrat(ptMaker(p1val,q2,p2))) - h);

%We choose initial conditions which have the correct energy one iterate
%forward
startPts = zeros(getnExtended(c),pointDensity);
parfor i = 1:pointDensity
    
    disp(i)
    %centerrads(i) = fzero(@(centerrad)refineHandle(centerrad,p1vals(i)),...
    %                      centerradguess)
    
    xopt = fsolve(@(x)multiRefineHandle(p1vals(i),x(1),x(2),h,...
              ptMaker,H2energy,revquadrat),[q2guess p2guess].',...
              optimset('TolX',1e-7,'TolFun',1e-7,'MaxFunEvals',2000));
           
%     xopt = fminunc(@(x)multiRefineHandle(p1vals(i),x(1),x(2),h,...
%                ptMaker,H2energy,revquadrat),[q2guess p2guess].',...
%                optimset('Display','iter'));       
    
    startPts(:,i) = ptMaker(p1vals(i),xopt(1),xopt(2));
end

xquad(:,1,:) = startPts;

%We map forward some iterates
for i = 2:iterates
    for j = 1:size(xquad,3)
        xquad(:,i,j) = revquadrat(squeeze(xquad(:,i-1,j)));
    end
end

manifoldPts = reshape(xquad(:,2:end,:),4,[]);

cplot(manifoldPts,c,'.k');

% %Now, we integrate. For speed, we enable caching.
% c = startCaching(c);
% 
% for i = 1:size(nontransitics,2)
%     [~,c] = integplot(linspace(0,timescale,numPts),...
%                       nontransitics(:,i),c,'bo');
% end
% for i = 1:size(transitics,2)
%     [~,c] = integplot(linspace(0,timescale,numPts),...
%                       transitics(:,i),c,'ro');
% end
% % for i = 1:size(manifold,2)
% %     [~,c] = integplot(linspace(0,timescale,numPts),...
% %                       manifold(:,i),c,'mo');
% % end
% 
% c = stopCaching(c);

c = cs(c,'s.o.v.dmode','position');

%We have to reset the coordinate system for now.
c = coordreset(c);

%Now, we integrate. For speed, we enable caching.
c = startCaching(c);

timescale = 2 * getSolarPeriod(c);
numPts = 1000;

nontransitics = cg(c,'lm.nontransit');
transitics = cg(c,'lm.transit');
manifold = cg(c,'lm.manifold');

integfig = figure
hold on
for i = 1:size(nontransitics,2)
    [p,c] = integplot(linspace(0,timescale,numPts),...
                      nontransitics(:,i),c,'b');
    p.Color(4) = 0.75;
    [p,c] = integplot(linspace(0,-timescale,numPts),...
                      nontransitics(:,i),c,'b');
    p.Color(4) = 0.75;
end
for i = 1:size(transitics,2)
    [p,c] = integplot(linspace(0,timescale,numPts),...
                      transitics(:,i),c,'r');
    p.Color(4) = 0.75;
    [p,c] = integplot(linspace(0,-timescale,numPts),...
                      transitics(:,i),c,'r');
    p.Color(4) = 0.75;
end
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

%Now, we plot the trajectory energies. We have to move back to the
%symplectic eigenbasis to do so, which in turn requires restarting the
%cache.

c = stopCaching(c);

% c = coordset(c,eigenbasis);
% 
% c = startCaching(c);
% 
% nontransitics = cg(c,'lm.nontransit');
% transitics = cg(c,'lm.transit');
% manifold = cg(c,'lm.manifold');
% 
% periods = 1;
% timescale = periods * getSolarPeriod(c);
% numPts = periods + 1;
% 
% figure
% hold on
% for i = 1:size(nontransitics,2)
%     [~,c] = energyplot(linspace(0,timescale,numPts),...
%                       nontransitics(:,i),c,H2energy,false,'c');
%     [~,c] = energyplot(linspace(0,-timescale,numPts),...
%                       nontransitics(:,i),c,H2energy,false,'c');              
% end
% for i = 1:size(transitics,2)
%     [~,c] = energyplot(linspace(0,timescale,numPts),...
%                       transitics(:,i),c,H2energy,false,'g');
%     [~,c] = energyplot(linspace(0,-timescale,numPts),...
%                       transitics(:,i),c,H2energy,false,'g');
% end
% % for i = 1:size(manifold,2)
% %     [~,c] = energyplot(linspace(0,timescale,numPts),...
% %                       manifold(:,i),c,H2energy,false,'m');
% %     [~,c] = energyplot(linspace(0,-timescale,numPts),...
% %                       manifold(:,i),c,H2energy,false,'m');
% % end
% 
% xticks(linspace(-timescale,timescale,2*numPts-1));
% 
% disp('CONDITION 3:')
% disp('Assess to see whether it holds visually.')
% disp('------------')
% 
% c = stopCaching(c);

function y = smootherabs(x,a)
    %no, this function isn't where you shave all the hair off your
    %abdominal muscles. It's a function similar to the absolute value
    %function, but smooth enough to be better behaved numerically. Found on 
    %https://math.stackexchange.com/questions/1284946/soft-absolute-value.
    
    y = a * abs(x)^3/(1+a*x^2);
end