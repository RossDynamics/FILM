%Obtains and performs some analysis on the monodromy matrix for the
%Earth-Moon BCM.

%c = Earth_Moon_Sun_BCM_Context;

c = coordreset(c);

disp('Switching to momentum coordinates...')
c = useMomentum(c);

%We now obtain the state transition matrices for each "piece"
c = cs(c,'p.thetam00',0);
stms = stm([0 getSolarPeriod(c)/4],cg(c,'lm.y0.t0'),c);
M1 = stms(:,:,2)

c = cs(c,'p.thetam00',pi/2);
stms = stm([0 getSolarPeriod(c)/4],cg(c,'lm.y0.t1'),c);
M2 = stms(:,:,2)

c = cs(c,'p.thetam00',pi);
stms = stm([0 getSolarPeriod(c)/4],cg(c,'lm.y0.t2'),c);
M3 = stms(:,:,2)

c = cs(c,'p.thetam00',3*pi/2);
stms = stm([0 getSolarPeriod(c)/4],cg(c,'lm.y0.t3'),c);
M4 = stms(:,:,2)

Mpieced = M4*M3*M2*M1

c = cs(c,'p.thetam00',0);
stms = stm([0 -getSolarPeriod(c)],cg(c,'lm.y0.t0'),c);
M = stms(:,:,2)

disp('M Entry Discrepancy: ')
disp(Mpieced - M)

peigenvals = eig(Mpieced,'vector')
eigenvals = eig(M,'vector')

disp('Reciprocal (Mpieced; Real eigenvalue > 1)')
disp(1/peigenvals(1))
disp('Error:')
disp(abs(peigenvals(2) - 1/peigenvals(1)))

disp('Reciprocal (M; Real eigenvalue > 1)')
disp(1/eigenvals(1))
disp('Error:')
disp(abs(eigenvals(2) - 1/eigenvals(1)))


%c = cs(c,'lm.monodromy',stms(:,:,2),2);

%disp('Monodromy matrix M:')
%disp(cg(c,'lm.monodromy'))