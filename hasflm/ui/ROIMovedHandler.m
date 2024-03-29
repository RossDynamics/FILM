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
extendedy0 = cg(c,'hasflm.s.extendedy0');
%We intentionally get the unextended n.
n = cg(c,'d.n');

%We get y0 from all of the roi's on the axes and from any extension values
%passed through
y0 = getIC(roiSet,extendedy0);

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
    else
        rethrow(exception)
    end
end

%We now plot on each axis.
modes = ["1-4" "2-3"];
for i = 1:2
    %We make sure to change the current axes and dimension mode
    c = cs(c,'s.o.v.dmode',modes(i));
    
    currentAxis = axisSet{i};
    axes(currentAxis);
    
    %We remove any objects with the deletion tag
    delete(findobj(currentAxis,'Tag','delete'));
    
    p = cplot(y,c,'b');
    p.Tag = 'delete';
    
    cplot(y0,c,'ok');
    
end

%We get the tolerance by calculating the distance between the start and end
%points on the trajectory, not counting extended phase space points
tolLabel.String = norm(y(1:n,1) - y(1:n,end));

end

