%Little bits of code for doing state transition tensor analysis that don't
%belong in their own scripts

%%
%Sets everything up
%mp.Digits(48);
load('C:\DATA\Virginia Tech\Research\CMDS Ecosystem\FILM_info\cTrue0-3 (2020_11_21 15_07_56 UTC).mat')
c = coordreset(c);
c = useMomentum(c);
integplot(linspace(0,getSolarPeriod(c),1000),mp(cg(c,'lm.y0')),c)
[phi,y,c] = stm(linspace(0,getSolarPeriod(c),3),mp(cg(c,'lm.y0')),c,[],2)

syms q1 q2 p1 p2
vars = [q1 q2 p1 p2].';

%%
%Verifies preservation of symplecticity

%BE SURE TO TURN THIS BACK OFF. IT CAUSES CMDS TO GIVE INCORRECT RESULTS
%DURING NUMERICAL INTEGRATION WHEN IT'S ON
sympref('FloatingPointOutput',true);  

%For the linear map:
disp('Verifying symplecticity for the linear map (should be J):')
phi{1}(:,:,end)'*Jmatrix(4)*phi{1}(:,:,end)

%For the quadratic map:
%Because gradient doesn't work on vectors, we have to use arrayfun.
%This isn't the correct way to verify symplecticity for state transition
%tensors, although I'll grant this tries to show the *map* is symplectic.
disp('Verifying symplecticity for the quadratic map (should be J):')
quadJacobianCell = arrayfun(@(expr)gradient(expr,vars).',...
               quadratlooping(vars,phi),'UniformOutput',false);

quadJacobianSymMat = cell2sym(quadJacobianCell);

simplify(quadJacobianSymMat.'*Jmatrix(4)*quadJacobianSymMat)

sympref('FloatingPointOutput',false);

%%
%Verifies preservation of symplecticity of STM's.
n = getnExtended(c);

phi1 = phi{1}(:,:,end);
phi2 = phi{2}(:,:,:,end);

Jm = Jmatrix(n);
disp('Verifying symplecticity for the quadratic map (should be J):')

vereq2 = zeros(n,n,n,'mp');
for i = 1:n
    for j = 1:n
        for k = 1:n
           for l = 1:n
               for m = 1:n
                   vereq2(i,l,m) = vereq2(i,l,m) + phi2(j,i,m)*Jm(j,k)*phi1(k,l) + ...
                                   phi1(j,i)*Jm(j,k)*phi2(k,l,m);
                               
                   if i == 1 && l == 1 && m == 1
                       disp(vereq2(i,l,m))
                   end
               end
           end
        end
    end
end

vereq2

%mp.Digits(34);
% 
vereq2new = ttt(ttt(tensor(phi2),tensor(mp(Jmatrix(n))),1,1),tensor(phi1),3,1) +...
ttt(ttt(tensor(phi1),tensor(mp(Jmatrix(n))),1,1),tensor(phi2),2,1)
 
%%
%Analyzes canonical planes in eigenbasis. Runs analyzeMonodromy first.
%Overwrites existing phi matrix in standard coordinates!
analyzeMonodromy
c = coordset(c,eigenbasis)
[phi,y,c] = stm(linspace(0,getSolarPeriod(c),3),mp(cg(c,'lm.y0')),c,[],2)

defineMapHandles
c = cs(c,'lm.phi',phi);
exploreMapsEigenbasis(2,1,c)
exploreMapsEigenbasis(2,2,c)