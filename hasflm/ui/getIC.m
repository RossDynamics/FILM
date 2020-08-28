function y0 = getIC(roiSet)
%GETIC Gets the full initial condition y0 for a trajectory from a set of
%roi's.

nhalf = numel(roiSet);

%We first construct the initial condition from all of the separate axes
for i = 1:nhalf
    y0([i,nhalf+i],:) = roiSet{i}.Position.';
end

end

