%A script for exploring the behaviors of the monodromy matrix map and the
%quadratic monodromy tensor map.

clear xlin xquad
close all;

linear = @(x)phi{1}(:,:,end)*x;
quadrat = @(x)(phi{1}(:,:,end)*x +...
          double(1/2*ttv(ttv(tensor(phi{2}(:,:,:,end)),x,1),x,1)));

numIter = 2;
pointDensity = 30;

xlin(:,1,:) = 1e-11*threesphere(pointDensity);
xquad(:,1,:) = 1e-11*threesphere(pointDensity);

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