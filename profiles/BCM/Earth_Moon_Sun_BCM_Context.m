function c = Earth_Moon_Sun_BCM_Context(epsilonomega)
%EARTH_MOON_SUN_BCM_CONTEXT A version of BCM_Context pre-populated with
%parameters for the Sun-perturbed Earth-Moon BCM. Supply an epsilonomega to
%control omegam0. If omegam0 = 1, the real value will be used.

c = BCM_Context;

c = cs(c,'p.mu', 0.0121505816234336);
c = cs(c,'p.mu0',328900.5499999991152436);
c = cs(c,'p.a0',388.8111430233511214);
c = cs(c,'p.omegam0',0.9251959855182896/epsilonomega);
c = cs(c,'p.thetam00',0);

end
