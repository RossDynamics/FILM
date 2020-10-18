%Obtains and performs some analysis on the monodromy matrix for the
%Earth-Moon BCM.

%c = Earth_Moon_Sun_BCM_Context;

c = coordreset(c);

%We set the initial condition for the Lagrange manifold
%c = cs(c,'lm.y0',[0.836688080819354   0.000000000000000  -0.000000000000002   0.000430728236819].');
disp('Lagrange manifold initial condition:')
disp(cg(c,'lm.y0'))

disp('Switching to momentum coordinates...')
c = useMomentum(c);

%We now obtain the forward and backward monodromy matrices
stms = stm([0 getSolarPeriod(c)],cg(c,'lm.y0'),c);
Mfd = stms(:,:,2);
stms = stm([0 -getSolarPeriod(c)],cg(c,'lm.y0'),c);
Mbd = stms(:,:,2);

disp('Forward monodromy matrix M:')
disp(Mfd)
disp('Backward monodromy matrix M:')
disp(Mbd)

%We now get the raw, unprocessed eigensystems
[fdeigenvecs,fdeigenvals] = eig(Mfd,'vector');
[bdeigenvecs,bdeigenvals] = eig(Mbd,'vector');

disp('Eigenvalues of the forward monodromy matrix:')
disp(fdeigenvals)

disp('Eigenvalues of the backward monodromy matrix:')
disp(bdeigenvals)

disp('Eigenvectors of the forward monodromy matrix:')
disp(fdeigenvecs)

disp('Eigenvectors of the backward monodromy matrix:')
disp(bdeigenvecs)

Cfd = getSymplecticBasis(Mfd);
disp('Forward symplectic basis Cfd:')
disp(Cfd)

Cbd = getSymplecticBasis(Mbd);
disp('Backward symplectic basis Cbd:')
disp(Cbd)

%%

%We demonstrate that the new eigenspace will be symplectic.
disp('C'' * J * C =')
disp(C' * Jmatrix(cg(c,'d.n')) * C)

disp('Switching origin to the Lagrange manifold initial condition...')
c = cs(c,'ac.origin',cg(c,'lm.y0'));

disp('Switching active basis to C...')
c = cs(c,'ac.basis',cg(c,'lm.C'));

coordget(c)

eigenbasis = coordget(c);

%We store the quadratic Hamiltonian in the context object
syms x(t) y(t) p1(t) p2(t) a sigma T
tillambda = log(sigma)/T;
theta = acos(a);
tilomega = theta / T;
H2 = tillambda * x * p1 + 1/2*tilomega*(y^2 + p2^2);

c = cs(c,'lm.H2',H2);

disp('Quadratic Hamiltonian: ')
disp(cg(c,'lm.H2'));

%Now, we set the parameters that are needed for the paramScan.
c = cs(c,'p.sigma',max(eigenvals));

%This expression finds the real parts of the complex eigenvalues. We will
%pick the first real part we find (since they'll be the same for both 
%eigenvalues)
a = real(eigenvals(~(eigenvals == real(eigenvals))));
c = cs(c,'p.a',a(1));

c = cs(c,'p.T',getSolarPeriod(c));

H2energy = getEnergyHandle(c,cg(c,'lm.H2'));