function exploreMaps(numIter,c)
%A script for exploring the behaviors of the monodromy matrix map and the
%quadratic monodromy tensor map. 10 is a good value for numIter at the
%CR3BP Lagrange point, whereas 2 is a good value for numIter at an L_1
%Lagrange manifold. c is a context object.

clear xlin xquad
%close all;

%Expects phi to be stored at the following location.
phi = cg(c,'lm.phi');

linear = @(x)phi{1}(:,:,end)*x;
quadrat = @(x)(phi{1}(:,:,end)*x +...
    double(1/2*ttv(ttv(tensor(phi{2}(:,:,:,end)),x,2),x,2)));

%numIter = 10;
pointDensity = 50;

xlin(:,1,:) = 1e-12*threesphere(pointDensity);
xquad(:,1,:) = 1e-12*threesphere(pointDensity);

for i = 2:numIter
    xlin(:,i,:) = linear(squeeze(xlin(:,i-1,:)));
    for j = 1:size(xquad,3)
        xquad(:,i,j) = quadrat(squeeze(xquad(:,i-1,j)));
    end
end

xlinr = reshape(xlin,4,[]);
xquadr = reshape(xquad,4,[]);

cplot(xlinr,c,'.b')
hold on
cplot(xquadr,c,'.r')

end