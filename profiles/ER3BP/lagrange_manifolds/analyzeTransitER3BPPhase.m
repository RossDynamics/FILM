%Used for analyzing transit using the symplectic eigenbasis.
%Run analyzeMonodromyER3BP first; call with analyzePhases.

%Let's set the eigenbasis first
c = coordset(c,eigenbasis);

%We get the start time
st = cg(c,'temp.startTime');

%To consider initial conditions in the upper halfplane, set halfplane equal
%to 1. To consider initial conditions in the lower halfplane, set
%halfplane equal to -1.
halfplane = 1;

%We create an array of initial conditions from within the
%symplectic eigenbasis along the bounding line. We set const (where p1 =
%const + q1) small and h = H2(x) small.

const = halfplane * 4e-5;
h = 1e-8;

%We also set a beginning for the array's q1's. This beginning coincides 
%with the p1 + q1 = 0 line, so we have q1 = -const / 2.

q1 = -const / 2;

%The other end of
%the array will be found dynamically; it will occur when the initial
%conditions become complex. We specify a small arrayStep for creating 
%subsequent q2's.

arrayStep = halfplane * 1e-6;

%We get the current initial condition for the LM. To be sure that a bug
%isn't occurring when accessing the phase, we access the raw values
lm_y0_str = cg(c,'temp.lm.y0',false);
lm_y0 = lm_y0_str.value

disp('Phase:')
lm_y0(5)

%Now, we build the array using ics_energy_boundary:
arrayics = ics_energy_boundaryER3BP(q1,const,h,...
                               cg(c,'p.sigma'),cg(c,'p.a'),cg(c,'p.T'),...
                               lm_y0(5));
                          
while true
    q1 = q1 + arrayStep;
    ic = ics_energy_boundaryER3BP(q1,const,h,...
                               cg(c,'p.sigma'),cg(c,'p.a'),cg(c,'p.T'),...
                               lm_y0(5));
    %If we get complex values, we know we've reached the boundary
    if ~isreal(ic)
        break;
    end
    arrayics = [arrayics ic];
end

disp('------------')
disp('CONDITION 2:')
disp('Last element of the initial conditions array''s q1*p1')
disp(arrayics(1,end)*arrayics(3,end));

disp('h/lambdatilde:')
disp(h / (1/2*pi*log(cg(c,'p.sigma'))));
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

% disp('Running analyzeParabolaER3BP...')
% analyzeParabolaER3BP
% disp('analyzeParabolaER3BP complete.')

lambdatil = log(cg(c,'p.sigma'))/(2*pi);
fplot(@(q_1)h/(lambdatil*q_1))

c = cs(c,'s.o.v.dmode','1');

hold on
limScale = 0.1;
line([0 0], [-limScale limScale]);
line([-limScale limScale], [0 0]);
cplot(cg(c,'lm.nontransit'),c,'.b');
cplot(cg(c,'lm.transit'),c,'.r');
%cplot(cg(c,'lm.manifold'),c,'or');

set(gca,'FontSize',20)
xlabel('$q_1$','interpreter','latex')
ylabel('$p_1$','interpreter','latex')
axis([-0.1 0.25 -0.05 0.3]*4e-4)
axis square

periods = 1;
timescale = periods * 2*pi;
numPts = 5;

nontransitics = cg(c,'lm.nontransit');
transitics = cg(c,'lm.transit');
manifold = cg(c,'lm.manifold');

%Now, we integrate. For speed, we enable caching.
%c = startCaching(c);

% for i = 1:size(nontransitics,2)
%     disp(i)
%     [~,c] = integplot(linspace(st,st+timescale,numPts),...
%                       nontransitics(:,i),c,'co');
% end
% for i = 1:size(transitics,2)
%     disp(i)
%     [~,c] = integplot(linspace(st,st+timescale,numPts),...
%                       transitics(:,i),c,'go');
% end
% for i = 1:size(manifold,2)
%     disp(i)
%     [~,c] = integplot(linspace(st,st+timescale,numPts),...
%                       manifold(:,i),c,'ro');
% end
%
%c = stopCaching(c);

c = cs(c,'s.o.v.dmode','position');

%We have to reset the coordinate system for now.
c = coordreset(c);

%Now, we integrate. For speed, we enable caching.
c = startCaching(c);

timescale = 1 * 2*pi;
numPts = 1000;

nontransitics = cg(c,'lm.nontransit')
transitics = cg(c,'lm.transit')
manifold = cg(c,'lm.manifold')

integfig = figure
hold on

[p,c] = integplot(linspace(st,st+timescale,numPts),...
                      cg(c,'temp.lm.y0'),c,'k');

for i = 1:size(nontransitics,2)
    [p,c] = integplot(linspace(st,st+timescale,numPts),...
                      nontransitics(:,i),c,'b');
    p.Color(4) = 0.75;
    [p,c] = integplot(linspace(st,st-timescale,numPts),...
                      nontransitics(:,i),c,'b');
    p.Color(4) = 0.75;
end
for i = 1:size(transitics,2)
    [p,c] = integplot(linspace(st,st+timescale,numPts),...
                      transitics(:,i),c,'r');
    p.Color(4) = 0.75;
    [p,c] = integplot(linspace(st,st-timescale,numPts),...
                      transitics(:,i),c,'r');
    p.Color(4) = 0.75;
 end
% for i = 1:size(manifold,2)
%     [p,c] = integplot(linspace(0,st+timescale,numPts),...
%                       manifold(:,i),c,'r');
%     p.Color(4) = 0.5;              
%     [p,c] = integplot(linspace(0,st-timescale,numPts),...
%                       manifold(:,i),c,'r');
%     p.Color(4) = 0.5;
% end

set(gca,'FontSize',20)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
axis equal

disp('CONDITION 1:')
disp('Assess to see whether it holds visually.')
disp('------------')

c = stopCaching(c);

drawnow

% %%
% %Now, we plot the trajectory energies. We have to move back to the
% %symplectic eigenbasis to do so, which in turn requires restarting the
% %cache.
% 
% c = coordset(c,eigenbasis);
% 
% c = startCaching(c);
% 
% nontransitics = cg(c,'lm.nontransit');
% transitics = cg(c,'lm.transit');
% manifold = cg(c,'lm.manifold');
% 
% periods = 1;
% timescale = periods * 2*pi;
% numPts = 2;
% 
% figure
% hold on
% for i = 1:size(nontransitics,2)
%     [~,c] = energyplot(linspace(st,st+timescale,numPts),...
%                       nontransitics(:,i),c,H2energy,false,'c');
%     [~,c] = energyplot(linspace(st,st-timescale,numPts),...
%                       nontransitics(:,i),c,H2energy,false,'c');              
% end
% for i = 1:size(transitics,2)
%     [~,c] = energyplot(linspace(st,st+timescale,numPts),...
%                       transitics(:,i),c,H2energy,false,'g');
%     [~,c] = energyplot(linspace(st,st-timescale,numPts),...
%                       transitics(:,i),c,H2energy,false,'g');
% end
% for i = 1:size(manifold,2)
%     [~,c] = energyplot(linspace(st,st+timescale,numPts),...
%                       manifold(:,i),c,H2energy,false,'r');
%     [~,c] = energyplot(linspace(st,st-timescale,numPts),...
%                       manifold(:,i),c,H2energy,false,'r');
% end
% 
% xticks(linspace(st-timescale,st+timescale,2*numPts-1));
% 
% disp('CONDITION 3:')
% disp('Assess to see whether it holds visually.')
% disp('------------')
% 
% c = stopCaching(c);