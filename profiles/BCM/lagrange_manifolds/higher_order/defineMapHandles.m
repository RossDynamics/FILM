%A short script that defines the function handles needed for much of this
%analysis. Be sure to have defined phi (the monodromy tensors of the
%Lagrange manifold) before running this code.

%Linear and quadratic maps (quadrat uses non-looping code. Doing so may 
%make it faster but also renders it dependent on the tensor toolbox package
%(for compatibility with older versions of MATLAB in which the initial
%conditions that we need were refined, we do not use tensorprod) and unable
%to accept symbolic inputs.
linearMap = @(x)phi{1}(:,:,end)*x;
quadrat = @(x)(phi{1}(:,:,end)*x +...
    double(1/2*ttv(ttv(tensor(phi{2}(:,:,:,end)),x,2),x,2)));