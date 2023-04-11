function c = Earth_Moon_Sun_CR3BP_Context(epsilonomega)
%EARTH_MOON_SUN_BCM_CONTEXT A version of BCM_Context pre-populated with
%parameters for the Sun-perturbed Earth-Moon CR3BP.
c = CR3BP_Context;

c = cs(c,'p.mu', 0.0121505816234336);

end
