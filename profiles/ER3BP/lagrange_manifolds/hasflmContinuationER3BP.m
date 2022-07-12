%Uses HASFLM and numerical continuation to attempt to find the 1e-0
%Lagrange manifold in the ER3BP.

%Let's create a context for the known Lagrange manifold.
%disp('Creating context with epsilonomega = 0...');
%cSet{1} = Earth_Moon_ER3BP_Context(0);

for i = 3:9
    
    disp('Iteration:')
    disp(i)
    
    c = cSet{i};
    
    if i > 1
        %Let's now run analyzeMonodromy.
        %disp('Running analyzeMonodromy...')
        %analyzeMonodromyER3BP
        cSet{i} = cs(cSet{i},'ac.origin',cg(cSet{i},'lm.y0.value',false));
        eigenbasis = coordget(cSet{i});
        cSet{i} = coordset(cSet{i},coordget(cSet{i-1}));
    else
        %Our initial guess is the L1 point location.
        xloc = 0.836915145386501754071018087383;
        cSet{i} = cs(cSet{i},'ac.origin',[xloc 0 0 0 0].');
        eigenbasis = coordget(cSet{i});
    end
    
    epsilone = (i+1)*0.1;

    %We now create a second context for a new epsilone.
    fprintf('Creating new context at %f\n',epsilone);
    cSet{i+1} = Earth_Moon_ER3BP_Context(epsilone);

    %if i > 1
    %We switch the old context's basis just in case analyzeMonodromy ever
    %changes (and doesn't leave the new context in the eigenbasis)
    disp('Switching contexts to the new Lagrange manifold eigenbasis...');

    cSet{i+1} = coordset(cSet{i+1},eigenbasis);
    %end
    
    clear HASFLM_y0;
    
    disp('Launching HASFLM...');
    
    %We insert a singularity detection event function to stop trajectories
    %before they get "stuck" next to a singularity.
    mu = cg(cSet{i+1},'p.mu');
    cSet{i+1} = cs(cSet{i+1},'s.i.odeopts',odeset('Events',...
                        @(t,y)singularityEventER3BP(t,y,mu,eigenbasis,cg(cSet{i+1},'p.e'))));
    
    cSet{i+1} = hasflm(2*pi,cSet{i+1},0);

    while ~exist('HASFLM_y0')
        pause;
    end
    
    disp('Setting the new guess...')
    cSet{i+1} = cs(cSet{i+1},'lm.y0',HASFLM_y0);
    
    %We remove the singularity detector after HASFLM finishes.
    cSet{i+1} = cs(cSet{i+1},'s.i.odeopts',odeset('Events',[]));

end