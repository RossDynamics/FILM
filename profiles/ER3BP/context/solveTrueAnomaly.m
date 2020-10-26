function c = solveTrueAnomaly(density,c)
%SOLVETRUEANOMALY Solves the true anomaly differential equation at density
%points and puts the resulting solution structure in c. Use getTrueAnomaly
%to access the true anomaly at a particular time. p.e (the eccentricity)
%and p.f0 (the true anomaly at t = 0) must be set first.

%The true anomaly over all time appears to be equal to the true anomaly
%from 0 to 2*pi times an offset (need to rigorously prove this!)
sol = integ(linspace(0,2*pi,density),cg(c,'p.f0'),c,...
        @(t,f)trueAnomalyFun(t,f,cg(c,'p.e')));

c = cs(c,'p.fsol',sol);

end

function fdot = trueAnomalyFun(~,f,e)
    fdot = ((1 + e*cos(f))^2)/((1-e^2)^(3/2));
end
