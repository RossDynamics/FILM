function c = BCM_Context()
%BCMCONTEXT A constructor/profile for the planar bicircular model that
%creates a context c.
c = Context(4);

syms x(t) y(t) vx(t) vy(t) p1(t) p2(t) t mu mu0 a0 omegam0 thetam00
c = variableNames([x y],[vx vy],[p1 p2],t,c);

%We input symbolic variables into the kinetic and potential energy
%functions.
T = BCM_T(x, y, vx, vy);
U = BCM_U(x, y, t, mu, mu0, a0, omegam0, thetam00);

c = cs(c,'d.T',T); 
c = cs(c,'d.U',U);

c = solveDynamics(c);

end

