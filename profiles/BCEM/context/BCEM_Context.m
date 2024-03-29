function c = BCEM_Context()
%ER3BP_CONTEXT A constructor/profile for the eccentric three-body problem
%that creates a context c. Run populateER3BPContext to set up the context's
%parameters.
c = Context(4);

syms x(t) y(t) vx(t) vy(t) p1(t) p2(t) t mu f(t) e mu0 a0 omegam0 thetam00
c = variableNames([x y],[vx vy],[p1 p2],t,c);

%We input symbolic variables into the kinetic and potential energy
%functions.
T = BCEM_T(x, y, vx, vy);
U = BCEM_U(x, y, t, mu, f, e, mu0, a0, omegam0, thetam00);

c = cs(c,'d.T',T); 
c = cs(c,'d.U',U);

c = solveDynamics(c);

end

