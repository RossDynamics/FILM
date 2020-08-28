function c = Earth_Moon_Sun_BCM_Context(epsilonomega)
%EARTH_MOON_SUN_BCM_CONTEXT A version of BCM_Context pre-populated with
%parameters for the Sun-perturbed Earth-Moon BCM. Supply an epsilonomega to
%control omegam0. If omegam0 = 1, the real value will be used.

c = BCM_Context;

c = cs(c,'p.mu',0.01215);
c = cs(c,'p.mu0',328900.54);
c = cs(c,'p.a0',388.81114);
c = cs(c,'p.omegam0',0.925195985520347/epsilonomega);
c = cs(c,'p.thetam00',0);

end
