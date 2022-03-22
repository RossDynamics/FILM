%Used for analyzing transit using the symplectic eigenbasis.
%Run analyzeMonodromy first.

%Let's set the eigenbasis first
c = coordset(c,eigenbasis);

%To consider initial conditions in the upper halfplane, set halfplane equal
%to 1. To consider initial conditions in the lower halfplane, set
%halfplane equal to -1.
halfplane = 1;

%We create an array of initial conditions from within the
%symplectic eigenbasis along the bounding line. We set const (where p1 =
%const + q1) small and h = H2(x) small.

const = halfplane * 5e-8;
h = 1e-10;

%We also set a beginning for the array's q1's. This beginning coincides 
%with the p1 + q1 = 0 line, so we have q1 = -const / 2.

q1 = -const / 2;

%The other end of
%the array will be found dynamically; it will occur when the initial
%conditions become complex. We specify a small arrayStep for creating 
%subsequent q2's.

arrayStep = halfplane * 2e-6;

%Now, we build the array using ics_energy_boundary:
arrayics = ics_energy_boundary(q1,const,h,...
                               cg(c,'p.sigma'),cg(c,'p.a'),cg(c,'p.T'));
                          
while true
    q1 = q1 + arrayStep;
    ic = ics_energy_boundary(q1,const,h,...
                               cg(c,'p.sigma'),cg(c,'p.a'),cg(c,'p.T'));
    %If we get complex values, we know we've reached the boundary
    if ~isreal(ic)
        break;
    end
    arrayics = [arrayics ic];
end

size(arrayics) 

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

timescale = 0.5 * getSolarPeriod(c);
numPts = 500;

nontransitics = cg(c,'lm.nontransit');
transitics = cg(c,'lm.transit');
manifold = cg(c,'lm.manifold');

integfig = figure
hold on

[lineSpec{1:size(nontransitics,2)}] = deal('b');
[lineSpec{size(nontransitics,2)+1:size(nontransitics,2)+size(transitics,2)}] = deal('r');
[lineSpec{size(nontransitics,2)+size(transitics,2)+1}] = deal('k');


integset = [nontransitics transitics cg(c,'lm.y0')];

% for i = 1:size(integset,2)
%     [sol,c] = integ(linspace(0,)
%     integset(:,i) =
% end

% set(gca,'FontSize',20)
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% axis equal
% 
% [F,c] = integanimate(linspace(0,timescale,numPts),...
%     integset,c,'bcmtransitfd',lineSpec);
% clf

set(gca,'FontSize',20)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
axis equal

[F,c] = integanimate(linspace(0,timescale,numPts),...
    integset,c,'bcmtransit',lineSpec);



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