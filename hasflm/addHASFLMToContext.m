function c = addHASFLMToContext(T,c,varargin)
%ADDHASFLMCONTEXT Sets up the HASFLM namespace of the context c with
%default properties and with period T. If an additional argument is
%supplied, each initial condition will incorprorate the extended phase 
%space initial conditions provided.

%Creates the root namespace
c.hasflm = struct;

%Creates a UI namespace
c.hasflm.ui = struct;

%Creates a settings namespace
c.hasflm.s = struct;
%The period of the desired periodic orbit; how long to integrate.
c.hasflm.s.T = Property(T,0);
%The number of points to use when integrating trajectories
c.hasflm.s.nPoints = Property(1000,0);

%Extended phase space coordinates, if applicable
if nargin >= 3
    extendedy0 = varargin{1}; 
else
    extendedy0 = [];
end
c.hasflm.s.extendedy0 = Property(extendedy0,0);

end

