%Obtains and performs some analysis on the monodromy matrix for the
%Earth-Moon ER3BP at a specific phase.

c = coordreset(c);

%We set the initial condition for the Lagrange manifold
%c = cs(c,'lm.y0',[0.836688080819354   0.000000000000000  -0.000000000000002   0.000430728236819].');
disp('Lagrange manifold initial condition:')
disp(cg(c,'temp.newy0'))

disp('Switching to momentum coordinates...')
c = useMomentum(c);

disp('Obtaining the monodromy matrix...')
%We now obtain the monodromy matrix (the state transition matrix after one
%period)
stms = stm([newTime 2*pi+newTime],cg(c,'temp.newy0'),c);
M = stms(1:4,1:4,2);
c = cs(c,'lm.monodromy',M,2);
Mfull = stms(:,:,2);
c = cs(c,'lm.monodromyfull',Mfull,2);

disp('Monodromy matrix M:')
disp(cg(c,'lm.monodromy'))

%We see what M' * J * M equals. It should equal J in momentum coordinates.
M = cg(c,'lm.monodromy');
disp('M'' * J * M =')
getnExtended(c)
disp(M' * Jmatrix(cg(c,'d.n')) * M)

%We now get the raw, unprocessed eigensystem
[eigenvecs,eigenvals] = eig(M,'vector');

disp('Eigenvalues of the monodromy matrix:')
disp(eigenvals)

disp('Eigenvectors of the monodromy matrix:')
disp(eigenvecs)

C = getSymplecticBasis(M,true);
%Because C is a basis, we store it using transform type 0.
c = cs(c,'lm.C',C,0);
disp('Symplectic basis C:')
disp(C)

%We demonstrate that the new eigenspace will be symplectic.
disp('C'' * J * C =')
disp(C' * Jmatrix(cg(c,'d.n')) * C)

%We use transformType = 0 because we're in a different basis.
disp('Switching origin to the Lagrange manifold initial condition...')
c = cs(c,'ac.origin',newy0,0);

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

c = cs(c,'p.T',2*pi);

H2energy = getEnergyHandle(c,cg(c,'lm.H2'));