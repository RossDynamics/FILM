%Some preparatory analysis before constructing a Lagrange manifold finder.

%Let's create a context for the known Lagrange manifold.
disp('Creating context with epsilonomega = 1e-1 (LM known)...');
c = Earth_Moon_Sun_BCM_Context(1e-1);

c = cs(c,'lm.y0',[0.836688080819354   0.000000000000000  -0.000000000000002   0.000430728236819].');

%Let's now run analyzeMonodromy.
disp('Running analyzeMonodromy...')
analyzeMonodromy

%We now create a second context for a new epsilonomega.
disp('Creating context with epsilonomega = 2e-1 (LM unknown)...');
c1 = Earth_Moon_Sun_BCM_Context(2e-1);

%We switch the old context's basis just in case analyzeMonodromy ever
%changes (and doesn't leave the new context in the eigenbasis)
disp('Switching contexts to the new Lagrange manifold eigenbasis...');
c = coordset(c,eigenbasis);
c1 = coordset(c1,eigenbasis);

disp('Setting the new guess...')
c1 = cs(c1,'lm.y0guess',cg(c,'lm.y0'));

%%

disp('Plotting guess in the new context...')

c1 = startCaching(c1);

close all;

pos = figure;
hold on;
c1 = cs(c1,'s.o.v.dmode','position');
[~,c1] = integplot(linspace(0,getSolarPeriod(c1),1000),...
                   cg(c1,'lm.y0guess'),c1,'.g');
cplot(cg(c1,'lm.y0guess'),c1,'.b');

%vel = figure;
%hold on;
% c1 = cs(c1,'s.o.v.dmode','velocity');
% [~,c1] = integplot(linspace(0,getSolarPeriod(c1),2),...
%                    cg(c1,'lm.y0guess'),c1,'.g');
% cplot(cg(c1,'lm.y0guess'),c1,'.b');

disp('Constructing offset in the new context...')
%We construct a hypersphere of offsets in an attempt to find a new IC.

offsetSphere = [];
spacing = logspace(0.5,-4,10);
for i = 1:size(spacing)
    offsetSphere = [offsetSphere spacing(i)*threesphere(7)];  
end

offsetSphere = offsetSphere + cg(c1,'lm.y0guess');

disp('Integrating offset in the new context...')
% figure(pos);
% hold on;
%c1 = cs(c1,'s.o.v.dmode','position');
fds = zeros(size(offsetSphere));
bds = zeros(size(offsetSphere));
for i = 1:size(offsetSphere,2)
    i
    [sol,c1] = integ(linspace(0,getSolarPeriod(c1),2),...
                   offsetSphere(:,i),c1); 
    fds(:,i) = sol.y(:,end);
end
[~,closeIndex] = sort(vecnorm(fds - offsetSphere));
close = offsetSphere(:,closeIndex(1))
%cplot(offsetSphere,c1,'.b');
%axis equal;

% figure(vel);
% hold on;
% c1 = cs(c1,'s.o.v.dmode','velocity');
% for i = 1:size(offsetSphere,2)
%     [~,c1] = integplot(linspace(0,getSolarPeriod(c1),2),...
%                    offsetSphere(:,i),c1,'.g'); 
% end
% for i = 1:size(offsetSphere,2)
%     [~,c1] = integplot(linspace(0,-getSolarPeriod(c1),2),...
%                    offsetSphere(:,i),c1,'.r'); 
% end
% cplot(offsetSphere,c1,'.b');
% axis equal;

c1 = stopCaching(c1);