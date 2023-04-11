function T = CR3BP_T(x, y, vx, vy)
%BCM_T The kinetic energy in the CR3BP expressed in velocity coordinates.

T = 1/2*((vx-y)^2 + (vy+x)^2);

end

