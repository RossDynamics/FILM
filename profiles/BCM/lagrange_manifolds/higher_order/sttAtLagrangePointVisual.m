c = Earth_Moon_ER3BP_e0_Context;
sttAtLagrangePoint

c = cs(c,'lm.phi',phi);

numIter = 10;

%Run with numIter rather large (10 seems to be a good value)
exploreMaps(numIter,c)

%We want to approximate the continuous time unstable manifold to the
%Lagrange point because it is our hypothesis that the unstable manifold is 
%shadowed by the linear and (more accurately) quadratic maps in a
%neighborhood of the equilibrium.
%
%The unstable eigenvector is unstabEigvec obtained using the following
%equations:
mu = cg(c,'p.mu');
xe = lppos(1);
mubar = mu*abs(xe-1+mu)^(-3)+(1-mu)*abs(xe+mu)^(-3);
a = 2*mubar + 1;
b = mubar - 1;
[eigvecs,eigvals] = eig([0 0   1 0
               0 0   0 1
               a  0  0 2
               0 -b -2 0])
unstabEigvec = real(eigvecs(:,2))

hold on;
c = useVelocity(c);
c = cs(c,'ac.origin',lppos);
integplot(linspace(0,endTime*numIter,1000),[1e-12*unstabEigvec;0],c,'k')
integplot(linspace(0,endTime*numIter,1000),-[1e-12*unstabEigvec;0],c,'k')
c = coordreset(c);