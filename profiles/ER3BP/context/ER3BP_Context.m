function c = ER3BP_Context()
%ER3BP_CONTEXT A constructor/profile for the eccentric three-body problem
%that creates a context c. Run populateER3BPContext to set up the context's
%parameters.
c = Context(4);

syms x(t) y(t) vx(t) vy(t) p1(t) p2(t) t mu f(t) e
c = variableNames([x y],[vx vy],[p1 p2],t,c);

%We input symbolic variables into the kinetic and potential energy
%functions.
T = ER3BP_T(x, y, vx, vy);
U = ER3BP_U(x, y, t, mu, f, e);

c = cs(c,'d.T',T); 
c = cs(c,'d.U',U);

c = solveDynamics(c);

end

