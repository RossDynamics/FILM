%Little bits of code for doing state transition tensor analysis that don't
%belong in their own scripts

%%
%Sets everything up
load('C:\DATA\Virginia Tech\Research\CMDS Ecosystem\FILM_info\cTrue0-3 (2020_11_21 15_07_56 UTC).mat')
integplot(linspace(0,getSolarPeriod(c),1000),cg(c,'lm.y0'),c)
[phi,y,c] = stm(linspace(0,getSolarPeriod(c),3),cg(c,'lm.y0'),c,[],2)

syms q1 q2 p1 p2
vars = [q1 q2 p1 p2].';

%%
%Verifies preservation of symplecticity

%BE SURE TO TURN THIS BACK OFF. IT CAUSES CMDS TO GIVE INCORRECT RESULTS
%DURING NUMERICAL INTEGRATION WHEN IT'S ON
sympref('FloatingPointOutput',true);  

%For the linear map:
disp('Verifying symplecticity for the linear map (should be J):')
phi{1}(:,:,end)*Jmatrix(4)*phi{1}(:,:,end)'

%For the quadratic map:
%Because gradient doesn't work on vectors, we have to use arrayfun
disp('Verifying symplecticity for the quadratic map (should be J):')
quadJacobianCell = arrayfun(@(expr)gradient(expr,vars).',...
               quadratlooping(vars,phi),'UniformOutput',false);

quadJacobianSymMat = cell2sym(quadJacobianCell);

simplify(quadJacobianSymMat*Jmatrix(4)*quadJacobianSymMat.')

sympref('FloatingPointOutput',false);

%%



%%
%Analyzes canonical planes in eigenbasis. Runs analyzeMonodromy first.
%Overwrites existing phi matrix in standard coordinates!
analyzeMonodromy
c = coordset(c,eigenbasis)
[phi,y,c] = stm(linspace(0,getSolarPeriod(c),3),cg(c,'lm.y0'),c,[],2)

defineMapHandles
c = cs(c,'lm.phi',phi);
exploreMapsEigenbasis(2,1,c)
exploreMapsEigenbasis(2,2,c)