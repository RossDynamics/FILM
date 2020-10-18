function circle = onesphere(n)
%ONESPHERE Creates a unit one-sphere (a circle) with point density n.

angles = linspace(0, 2*pi, n+1);

angles = angles(1:(end-1));

circle = [cos(angles).' sin(angles).']

end

