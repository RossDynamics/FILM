function c = Earth_Moon_ER3BP_Context(varargin)
%EARTH_MOON_ER3BP_Context Creates an ER3BP context for the Earth-Moon
%system.

%We can optionally rescale the perturbation
if nargin > 0
    escale = varargin{1};
else
    escale = 1;
end

c = ER3BP_Context;

c = cs(c,'p.mu',0.0121505816234336);
c = cs(c,'p.e',0.0549006*escale);

syms f(t) e
fEqn = diff(f,t) == (1 + e*cos(f))^2/((1-e^2)^(3/2));

c = addIntegratedParameter(f,fEqn,c);

end

