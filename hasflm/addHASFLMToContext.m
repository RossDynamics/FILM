function c = addHASFLMToContext(T,c)
%ADDHASFLMCONTEXT Sets up the HASFLM namespace of the context c with
%default properties and with period T.

%Creates the root namespace
c.hasflm = struct;

%Creates a UI namespace
c.hasflm.ui = struct;

%Creates a settings namespace
c.hasflm.s = struct;
%The period of the desired periodic orbit; how long to integrate.
c.hasflm.s.T = Property(T,0);
%The number of points to use when integrating trajectories
c.hasflm.s.nPoints = Property(100,0);

end

