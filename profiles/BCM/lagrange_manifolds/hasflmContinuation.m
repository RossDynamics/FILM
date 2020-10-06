%Uses HASFLM and numerical continuation to attempt to find the 1e-0
%Lagrange manifold.

%Let's create a context for the known Lagrange manifold.
disp('Creating context with epsilonomega = 1e-1 (LM known)...');
cSet{1} = Earth_Moon_Sun_BCM_Context(1e-1);

cSet{1} = cs(cSet{1},'lm2.y0',[1.1 0 0 0].');

%cSet{8} = Earth_Moon_Sun_BCM_Context(8e-1);

%cSet{8} = cs(cSet{8},'lm.y0',[0.821153761674520  -0.000008162988573   0.000033569946083   0.143822404285661].');

%cSet{2} = Earth_Moon_Sun_BCM_Context(2e-1);

%cSet{2} = cs(cSet{2},'lm.y0',[0.836616983019289   0.000000005805273   0.000000001075985   0.000996805459492].');

for i = 1:9
    
    disp('Iteration:')
    disp(i)
    
    c = cSet{i};
    
    %Let's now run analyzeMonodromy.
    disp('Running analyzeMonodromy...')
    analyzeMonodromy
    
    epsilonomega = (i+1)*0.1;

    %We now create a second context for a new epsilonomega.
    fprintf('Creating new context at %f\n',epsilonomega);
    cSet{i+1} = Earth_Moon_Sun_BCM_Context(epsilonomega);

    %We switch the old context's basis just in case analyzeMonodromy ever
    %changes (and doesn't leave the new context in the eigenbasis)
    disp('Switching contexts to the new Lagrange manifold eigenbasis...');
    c = coordset(c,eigenbasis);
    cSet{i+1} = coordset(cSet{i+1},eigenbasis);
    
    clear HASFLM_y0;
    
    disp('Launching HASFLM...');
    
    %We insert a singularity detection event function to stop trajectories
    %before they get "stuck" next to a singularity.
    mu = cg(cSet{i+1},'p.mu');
    cSet{i+1} = cs(cSet{i+1},'s.i.odeopts',odeset('Events',...
                        @(t,y)singularityEvent(t,y,mu,eigenbasis)));
    
    cSet{i+1} = hasflm(getSolarPeriod(cSet{i+1}),cSet{i+1});

    while ~exist('HASFLM_y0')
        pause;
    end
    
    disp('Setting the new guess...')
    cSet{i+1} = cs(cSet{i+1},'lm.y0',HASFLM_y0);
    
    %We remove the singularity detector after HASFLM finishes.
    cSet{i+1} = cs(cSet{i+1},'s.i.odeopts',odeset('Events',[]));

end