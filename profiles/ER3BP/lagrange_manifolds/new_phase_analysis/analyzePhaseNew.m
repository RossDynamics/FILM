c = coordreset(c);

c = useMomentum(c);

%We calculate f from the new time while integrating to the new initial
%condition.
newTime = 0;

if newTime
    %We integrate the initial condition to the new time. Be sure not to go too
    %far forward unless you have an extremely accurate initial condition
    [sol,c] = integ([0 newTime],cg(c,'lm.y0'),c);

    %We get the new initial condition
    newy0 = sol.y(:,end);
    c = cs(c,'temp.newy0',newy0);

    disp('New initial condition:')
    disp(newy0)
else
    newy0 = cg(c,'lm.y0');
    c = cs(c,'temp.newy0',newy0);
end

analyzeMonodromyPhaseNew
analyzeTransitPhaseNew

%integplot(linspace(newTime, newTime + 2*pi,500),newy0,c)