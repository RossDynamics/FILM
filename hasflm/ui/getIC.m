function y0 = getIC(roiSet,extensions)
%GETIC Gets the full initial condition y0 for a trajectory from a set of
%roi's. If the phase space was extended, initial conditions for the
%extended values are also passed in.

nhalf = numel(roiSet);

%We first construct the initial condition from all of the separate axes
y0([1,4],:) = roiSet{1}.Position.';
y0([2,3],:) = roiSet{2}.Position.';

y0 = [y0; extensions];
end

