%Uses HASFLM and numerical continuation to attempt to find the 1e-0
%Lagrange manifold.

%Let's create a context for the known Lagrange manifold.
%disp('Creating context with epsilonomega = 1e-1 (LM known)...');
%cSet{1} = Earth_Moon_Sun_BCM_Context(1e-1);

%cSet{1} = cs(cSet{1},'lm.y0',[0.836688080819354   0.000000000000000  -0.000000000000002   0.000430728236819].');

cSet{4} = Earth_Moon_Sun_BCM_Context(4e-1);

cSet{4} = cs(cSet{4},'lm.y0',[0.836349913728578   0.000000037601371  -0.000000070846718   0.002605328602630].');

for i = 4:9
    
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
    
    cSet{i+1} = hasflm(getSolarPeriod(cSet{i+1}),cSet{i+1});

    while ~exist('HASFLM_y0')
        pause;
    end
    
    disp('Setting the new guess...')
    cSet{i+1} = cs(cSet{i+1},'lm.y0',HASFLM_y0);

end