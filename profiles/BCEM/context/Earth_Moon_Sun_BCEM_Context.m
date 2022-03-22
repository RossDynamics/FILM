function c = Earth_Moon_Sun_BCEM_Context()
%EARTH_MOON_BCEM_Context Creates an BCEM context for the Earth-Moon
%system. The BCEM is the BCM, but with the primaries moving as in the
%ER3BP.

c = BCEM_Context;

c = cs(c,'p.mu',0.0121505816234336);
c = cs(c,'p.e',0.0549006);
c = cs(c,'p.mu0',328900.5499999991152436);
c = cs(c,'p.a0',388.8111430233511214);
c = cs(c,'p.omegam0',0.9251959855182896);
c = cs(c,'p.thetam00',0);

syms f(t) e
fEqn = diff(f,t) == (1 + e*cos(f))^2/((1-e^2)^(3/2));

c = addIntegratedParameter(f,fEqn,c);

end

