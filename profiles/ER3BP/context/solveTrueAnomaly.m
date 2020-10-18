function c = solveTrueAnomaly(t1,c)
%SOLVETRUEANOMALY Solves the true anomaly differential equation and puts
%the resulting solution structure in c. Use getTrueAnomaly to access the
%true anomaly at a particular time. p.e (the eccentricity) and p.f0 (the
%true anomaly at t = 0) must be set first.

%The true anomaly
sol = integ([0 t1],cg(c,'p.f0'),c,@(t,f)trueAnomalyFun(t,f,cg(c,'p.e')));

%The true anomaly solution object is saved to p.fsol.
c = cs(c,'p.fsol',sol);

end

function fdot = trueAnomalyFun(~,f,e)
    fdot = ((1 + e*cos(f))^2)/((1-e^2)^(3/2));
end
