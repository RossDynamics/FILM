%Attempts to find the cubic Hamiltonian for the quadratic map built from
%the linear and quadratic monodromy tensors contained in phi.

syms q1 p1 q2 p2
%syms q p

% phi = {};
% phi{1}(:,:,1) = eye(2);
% phi{2}(:,:,1,1) = [40  10
%                 -6 -40];
% phi{2}(:,:,2,1) = [10  48
%                -40 -10];             
% 
% phi = {};
% phi{1}(:,:,1) = zeros(2);
% phi{2}(:,:,1,1) =  [0  0
%                    -1  0];
% phi{2}(:,:,2,1) = [0 1
%                    0 0];

defineMapHandles

HOrder = 2;
workingOrder = 3;
% HOrder = 3;
% workingOrder = 7;

%vars = [q p].';
vars = [q1 q2 p1 p2].';
HMonomials = getMonomials(vars,HOrder);
workingMonomials = getMonomials(vars,workingOrder);
%HMonomials = [getMonomials(vars,HOrder); [1/q 1/p q/p p/q 1/q^2 1/p^2 1].'];
%workingMonomials = [getMonomials(vars,workingOrder); [1/q 1/p q/p p/q 1/q^2 1/p^2 1].'];
map = linearMap
%map = @(x)quadratlooping(x,phi)
%sigma = 1e+12;
% map = @(x)[sigma*x(1)
%          cos(pi/6)*x(2)+sin(pi/6)*x(4)
%          (1/sigma)*x(3)
%         -sin(pi/6)*x(2)+cos(pi/6)*x(4)]

sympref('FloatingPointOutput',true);      
      
mapH = computeMapHamiltonian(map,HMonomials,workingMonomials,vars)

bnum = numel(children(mapH));
syms b [bnum 1]

B = eye(numel(bnum));

mapHex = subs(mapH,b,ones(bnum,1))
%mapHex = subs(mapH,b,B(:,1))
double(subs(mapHex,vars,[5 2 3 4].'))
double(subs(mapHex,vars,map([5 2 3 4].')))

sympref('FloatingPointOutput',false);