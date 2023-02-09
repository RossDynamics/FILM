%Calculates the linear and quadratic state transition tensors at the
%unperturbed Lagrange point in the ER3BP for a desired amount of time and
%compares them to a more analytical method of calculation as a check. Use
%velocity coordinates.

c = useVelocity(c);

%We get the location of the L1 point
lppos = [findlp(cg(c,'p.mu'),0.5);0;0;0;0];

%We calculate the STT's
endTime = 1;
[phi,y,c] = stm(linspace(0,endTime,3),lppos,c,[],2);
phi = {phi{1}(1:4,1:4,:) phi{2}(1:4,1:4,1:4,:)};

disp('Linear STT:')
disp(phi{1}(:,:,3))

%Now, we attempt to calculate via the analytical method. 
%We first get the equations of motion and function handles for calculating
%the Jacobians:
eqns = cg(c,'d.eqns');
J1 = getJacobianHandle(eqns,c,1);
J2 = getJacobianHandle(eqns,c,2);

%The right hand side of the equations of motion at the L1 point is equal to
%zero, so we can substitute the constant values in lppos into the Jacobian
%function J for computing the linear STT. Thus, the STM is of the form
%\dot{phi} = ... J1(lppos) * phi. But, given the initial condition phi(0) =
%eye(n), the solution is phi(t) = e^(J(lppos) * t):
phi1_recalc = expm(J1(0,[lppos;0])*endTime);
disp('Linear STT recalculated using a different method:')
disp(phi1_recalc)

%Now, we calculate the quadratic STT.
%At the Lagrange point, we notice something very interesting about the form
%of J2: the only nonzero elements (ignoring those associated with the
%phase space extension) are 311, 421, 412, and 322:
disp('Second Jacobian tensor:')
disp(J2(0,[lppos;0]))

%Similarly, we see that most of the elements of J1 are zero:
disp('First Jacobian tensor:')
disp(J1(0,[lppos;0]))

%We do a calculation by hand (see Research Notebook II) that permits the
%determination of the equations of motion. 

disp('Quadratic STT:')
disp(phi{2}(:,:,:,3))

%disp('Quadratic STT recalculated using a different method:')
%disp(phi2_recalc)

