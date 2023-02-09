%For this analysis, be sure that a context object with an appropiate
%manifold has already been loaded. Also, be sure to set the period of the
%Lagrange manifold in the variable period.

%period = 2*pi;
period = getSolarPeriod(c);

%We calculate the monodromy tensors
[phi,y,c] = stm(linspace(0,period,3),cg(c,'lm.y0'),c,[],2)

%We have to process phi to remove the extended phase space variables
phi{1} = phi{1}(1:4,1:4,:);
phi{2} = phi{2}(1:4,1:4,1:4,:);
c = cs(c,'lm.phi',phi);

linear = @(x)phi{1}(:,:,end)*x;
quadrat = @(x)(phi{1}(:,:,end)*x +...
    double(1/2*ttv(ttv(tensor(phi{2}(:,:,:,end)),x,2),x,2)));

numIter = 2;

%Run with numIter rather small (2 seems to be a good value)
exploreMaps(numIter,c)

%We get the unstable eigenvector from the monodromy matrix, which is just
%the eigenvector corresponding to the real eigenvalue > 1 from
%the linear STT:
[eigvecs,eigvals] = eig(phi{1}(:,:,end))
unstabEigvec = real(eigvecs(:,1))

%displacements = 1e-12*unstabEigvec;
displacements = logspace(-14,-8,100).*unstabEigvec;
if cg(c,'d.n') == 5
    displacements = [displacements; zeros(size(displacements,2),1)];
end

hold on;

c = useVelocity(c);
c = cs(c,'ac.origin',cg(c,'lm.y0'));

c = startCaching(c);

for i = 1:size(displacements,2)
%     cplot(linear(cg(c,'lm.y0')+displacements(:,i)),c,'b.')
%     cplot(linear(cg(c,'lm.y0')-displacements(:,i)),c,'b.')
%     cplot(quadrat(cg(c,'lm.y0')+displacements(:,i)),c,'r.')
%     cplot(quadrat(cg(c,'lm.y0')-displacements(:,i)),c,'r.')

    integplot(linspace(0,period*(numIter-1),numIter),cg(c,'lm.y0')+...
                       displacements(:,i),c,'ko')
    integplot(linspace(0,period*(numIter-1),numIter),cg(c,'lm.y0')-...
                       displacements(:,i),c,'ko')
end
              
integplot(linspace(0,period,1000),cg(c,'lm.y0'),c,'g')
axis equal

c = stopCaching(c);

c = coordreset(c);