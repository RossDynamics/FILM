function f = getTrueAnomaly(t,fvals,tvals)
%GETTRUEANOMALY Retrieves the current true anomaly f at time t from the
%solution fvals from 0 to 2*pi at times tvals.

%Because f seems to be able to be expressed by "gluing" f between 0 and
%2*pi together, we calculate f at an arbitrary time from f between 0 and
%2*pi.
tclamped = mod(t,2*pi);

%We find the closest clamped value
[~,i] = min(abs(tclamped - tvals));

f = fvals(:,i);

end

