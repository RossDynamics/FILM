function exploreMapsEigenbasis(numIter,canPlane,c)
%A script for exploring the behaviors of the monodromy matrix map and the
%quadratic monodromy tensor map in the linear eigenbasis. canPlane, the
%number of the canonical plane in which to generate coordinates, must be
%specified. 10 is a good value for numIter at the CR3BP Lagrange point,
%whereas 2 is a good value for numIter at an L_1 Lagrange manifold. c is a
%context object.

clear xlin xquad
%close all;

%Expects phi to be stored at the following location.
phi = cg(c,'lm.phi');

linear = @(x)phi{1}(:,:,end)*x;
quadrat = @(x)(phi{1}(:,:,end)*x +...
    double(1/2*ttv(ttv(tensor(phi{2}(:,:,:,end)),x,2),x,2)));

%Builds out coordinates that lie purely within the specified canonical
%planes
pointDensity = 50;
[q1Pts,p1Pts] = meshgrid(linspace(-1,1,pointDensity),...
                         linspace(-1,1,pointDensity));
q1p1Pts = [reshape(q1Pts,1,[]) 
           reshape(p1Pts,1,[])];
               
startPts = zeros(4,numel(q1Pts));
startPts(canPlane,:) = q1p1Pts(1,:);
%startPts(canPlane+cg(c,'d.n')/2,:) = q1p1Pts(2,:);

% [q1Pts,q2Pts,p2Pts] = meshgrid(linspace(-1,1,pointDensity),...
%                          linspace(-1,1,pointDensity),...
%                          linspace(-1,1,pointDensity));
% startPts = [reshape(q1Pts,1,[])
%             reshape(q2Pts,1,[])
%             zeros(1,numel(q1Pts))
%             reshape(p2Pts,1,[])];

epsilon = 1e-12;

xlin(:,1,:) = epsilon*startPts;
xquad(:,1,:) = epsilon*startPts;

for i = 2:numIter
    xlin(:,i,:) = linear(squeeze(xlin(:,i-1,:)));
    for j = 1:size(xquad,3)
        xquad(:,i,j) = quadrat(squeeze(xquad(:,i-1,j)));
    end
end

xlinr = reshape(xlin,4,[]);
xquadr = reshape(xquad,4,[]);

figure;

subplot(2,2,1);
c = cs(c,'s.o.v.dmode','1')

cplot(xlinr,c,'.b')
hold on
cplot(xquadr,c,'.r')

xlabel('$q_1$','Interpreter','Latex')
ylabel('$p_1$','Interpreter','Latex')

title('q1-p1 Canonical Plane')

subplot(2,2,2);
c = cs(c,'s.o.v.dmode','2')

cplot(xlinr,c,'.b')
hold on
cplot(xquadr,c,'.r')

xlabel('$q_2$','Interpreter','Latex')
ylabel('$p_2$','Interpreter','Latex')

c = cs(c,'s.o.v.dmode','position')

title('q2-p2 Canonical Plane')

subplot(2,2,3);
plot3(xquadr(1,:),xquadr(2,:),xquadr(3,:),'.')
hold on
plot3(xlinr(1,:),xlinr(2,:),xlinr(3,:),'.')

xlabel('$q_1$','Interpreter','Latex')
ylabel('$q_2$','Interpreter','Latex')
zlabel('$p_1$','Interpreter','Latex')

title('q_1-q_2-p_1 Space')

subplot(2,2,4);
plot3(xquadr(1,:),xquadr(3,:),xquadr(4,:),'.')
hold on
plot3(xlinr(1,:),xlinr(3,:),xlinr(4,:),'.')

xlabel('$q_1$','Interpreter','Latex')
ylabel('$p_1$','Interpreter','Latex')
zlabel('$p_2$','Interpreter','Latex')

title('q_1-p_1-p_2 Space')

sgtitle(['Mapping initial conditions lying purely within the q'...
       int2str(canPlane) '-p' int2str(canPlane) ' canonical plane'])
   

end