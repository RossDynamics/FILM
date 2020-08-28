function c = hasflm(T,c)
%HASFLM HASFLM is the Human-Assisted Solver For Lagrange Manifolds, a
%utility for gamifying the search for Lagrange manifolds.
%Pass in a context object c to use; also pass in a period T over which to
%integrate.

%Adds HASFLM properties to the current context.
c = addHASFLMToContext(T,c);

%We store the old cache state so that we can turn it back off if needed.
cacheWasOn = isCaching(c);
c = startCaching(c);

%We cache eqnsHandle manually because we can't easily retrieve it from the
%handler functions.
[~,c] = getFromCache(c,'ca.eqnsHandle',@defaultHandle,c);

[axisSet,pointSet] = createUI(c);

if ~cacheWasOn
    c = stopCaching(c);
end

end

