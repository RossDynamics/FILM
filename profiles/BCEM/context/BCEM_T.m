function T = BCEM_T(x, y, vx, vy)
%ER3BP_T The kinetic energy in the ER3BP expressed in velocity coordinates.

T = 1/2*((vx-y)^2 + (vy+x)^2);

end