function xcoord = findlp(mu,guess)
%FINDLP Computes the x coordinate of the collinear Lagrange point nearest
%to the initial guess guess
xcoord = fzero(@(x)(-x + (mu * (-1 + mu + x))/abs(-1 + mu + x)^3 - ...
                   ((-1 + mu)*(mu + x))/abs(mu + x)^3),guess);
end

