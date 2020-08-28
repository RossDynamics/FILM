function T = getSolarPeriod(c)
%GETSOLARPERIOD Calculates the solar period from the omegam0 parameter in
%the context c.

T = 2*pi/cg(c,'p.omegam0');

end

