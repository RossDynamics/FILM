function ROIMovedHandler(axisSet,roiSet,tolLabel,c)
%ROIMovedHandler Handles ROI events. src and evt should be accessible in
%function handles passed to addlistener. num identifies the qi - qdoti/pi
%coordinate plane for which the current handler is reponsible. axisSet
%contains all axes, and roiSet contains all roi points (including this
%one). c is a CMDS context object. All arguments except src and evt should
%be set before passing in the function handle.

nhalf = numel(roiSet);

T = cg(c,'hasflm.s.T');
nPoints = cg(c,'hasflm.s.nPoints');

%We get y0 from all of the roi's on the axes
y0 = getIC(roiSet);

tspan = linspace(0,T,nPoints);
sol = integ(tspan,y0,c);
try
    y = deval(sol,tspan);
catch exception
    %If an event function or something similar prematurely ends
    %integration, we have to get what we can
    if strcmp(exception.identifier,'MATLAB:deval:SolOutsideInterval')
        tspan = linspace(0,sol.x(end),nPoints);
        y = deval(sol,tspan);
        disp('Unable to integrate over the full timespan required.')
    end
end

%We now plot on each axis.
for i = 1:nhalf
    %We make sure to change the current axes and dimension mode
    c = cs(c,'s.o.v.dmode',int2str(i));
    
    currentAxis = axisSet{i};
    axes(currentAxis);
    
    %We remove any objects with the deletion tag
    delete(findobj(currentAxis,'Tag','delete'));
    
    p = cplot(y,c,'b');
    p.Tag = 'delete';
    
    cplot(y0,c,'ok');
end

%We get the tolerance by calculating the distance between the start and end
%points on the trajectory
tolLabel.String = norm(y(:,1) - y(:,end));


end

