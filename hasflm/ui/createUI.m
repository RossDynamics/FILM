function [axisSet,roiSet,tolLabel] = createUI(c,varargin)
%CREATEUI Creates a HASFLM UI for the context object c. Returns cell arrays
%containing the axes and roi points created. Also returns UI
%elements. Will take initial conditions for phase space extensions as an
%optional argument. This version of the UI has been rewritten to be more
%restrictive but efficient; it is designed only to find initial conditions
%of the form (q1,0,0,p2).

%We need to know the unextended phase space dimension so that we can know 
%how many axes to create
n = cg(c,'d.n');

close all;

tiledlayout('flow');

%We first create axes and ROI objects
axisSet{1} = nexttile;
    
title('$q_1, \dot{q}_2$ or $q_1, p_2$','Interpreter','latex')

%We obtain a drawpoint and add a listener
roiSet{1} = drawpoint('Position',[0 0]);

axis(1e-1*[-1 1 -1 1]);

hold on;

axisSet{2} = nexttile;
    
title('$q_2, \dot{q}_1$ or $q_2, p_1$','Interpreter','latex')

%We obtain a drawpoint and add a listener
roiSet{2} = drawpoint('Position',[0 0]);

axis(1e-1*[-1 1 -1 1]);

hold on;

%We now create a special axis, the statusAxis, for viewing data
statusAxis = nexttile([1 2]);
axis([-1 1 -1 1]);
text(-0.5,0,'Tolerance of current periodic orbit candidate:');
tolLabel = text(-0.5,-0.2,'Awaiting integration...');
set(gca,'visible','off')

tb = axtoolbar(statusAxis);
btn = axtoolbarbtn(tb,'push');
btn.Tooltip = 'Get Initial Condition';
if nargin >= 2
    extendedy0 = varargin{1}; 
else
    extendedy0 = [];
end
btn.ButtonPushedFcn = @(~,~)assignin('base','HASFLM_y0',getIC(roiSet,...
                                                        extendedy0));

%We create the listeners in a separate loop so that we can provide the
%entire axisSet, roiSet, and relevant UI elements to each listener.

for i = 1:n/2 
    listenerSet{i} = addlistener(roiSet{i},'ROIMoved',...
    @(~,~)ROIMovedHandler(axisSet,roiSet,tolLabel,c));
end

end

