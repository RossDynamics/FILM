function [mapH,mapHspace] = computeMapHamiltonian(map,HMonomials,...
                                            workingMonomials,vars,varargin)
% COMPUTEMAPHAMILTONIAN For the symplectic map provided, attempts to
% determine what subspace of the vector space of the monomials provided is
% conserved under the action of the map; that is, what Hamiltonians
% (conserved quantities) can be expressed as a linear combination of the
% provided monomials). map should be specified as a function handle and
% monomials should be specified as a symbolic column vector. The phase
% space variables used in monomials should also be provided as a symbolic
% column vector. For compatibility with CMDS, using the [q; p] sort order
% for expressions is highly recommended. mapH is the Hamiltonian; an extra
% output argument mapHspace allows direct access to the results of the
% nullspace computation. An optional input argument allows you to specify
% the magnitude threshold for terms; terms with coefficients smaller than
% this magnitude will be culled. This threshold defaults to 1e-40. A second
% optional input argument allows you to specify the tolerance within which
% singular values will be considered close enough to 0. This threshold
% defaults to 1e-3.

if nargin >= 5
    cullingThreshold = varargin{1};
else
    cullingThreshold = 1e-40;
end

if nargin >= 6
    svdThreshold = varargin{2};
else
    svdThreshold = 1e-3;
end

%We use the sorting orders specified in the inputs themselves, so the
%below lines of code are commented out.
%monomials = sort(monomials);
%vars = sort(vars);

%These are the coefficients out of which Hx will be formed. We form them
%initially as a column vector.
syms a [numel(HMonomials) 1]

%The following is the form that the Hamiltonian must have, as a linear
%combination of the coefficients and the monomials:
Hx = a.'*HMonomials;

%Now, we get the vector corresponding to an iterate of the map:
Lambdax = expand(map(vars));

%We create HLambdax, the Hamiltonian after one iterate, by substituting in
%Lambdax for vars. We expand to avoid problems and then collect terms so
%that everything is written in terms of the monomials:
HLambdax = collect(expand(subs(Hx,vars,Lambdax)),workingMonomials);

%Now, we break each term into a separate expression. May behave strangely
%if there aren't multiple terms. WARNING: This code assumes that children
%returns a vector, not a cell array, because it was written for R2020a.
%However, this behavior was changed in R2020b, and so it probably won't
%work in R2020b until I add a version case handler.
HLxterms = children(HLambdax);

A = sym([]);

%We are forced to loop in order to preserve the monomial order. I find this
%fact rather annoying, because with good symbolic sorting functionality 
%in MATLAB this code could probably be written without a loop.
parfor i = 1:numel(workingMonomials)
    
    a = sym('a',[numel(HMonomials) 1]);
    
    %We divide out the monomial from all the terms...
    HLxtermsdiv = simplify(HLxterms/workingMonomials(i));
    %...and then we identify the term in which none of the phase space
    %variables are present (that is, the term where there is only a
    %coefficient left).
    coe = HLxtermsdiv(~has(HLxtermsdiv,vars));
    
    if isempty(coe)
        coe = 0;
        A(i,:)=zeros(1,numel(HMonomials));
        continue;
    end
    
    %Now, we break apart the coefficient into its constituent parts,
    %but this time using the a(i) variables, only if there the coefficient
    %is a summation of more than one term
    if hasSymType(coe,'plus')
        coeterms = children(coe);
    else
        coeterms = coe;
    end
    Arow = [];
    for j = 1:numel(HMonomials)
        coetermsdiv = simplify(coeterms/a(j));
        
        coecoe = coetermsdiv(~has(coetermsdiv,a));
        
        if isempty(coecoe)
            coecoe = 0;
        end
        
        Arow = [Arow coecoe];
    end
  
    A(i,:) = Arow;
     
    disp(i)
end

%We now try to solve the equation HLambdax = Hx once it has been cast in
%this form
%mapHspace = null(A - eye(numel(monomials)));

%We define a matrix corresponding to the original value of Hx and then pad
%it
Hxeye = eye(numel(HMonomials));
if any(size(Hxeye) < size(A))
    Hxeye(size(A,1),size(A,2)) = 0;
end

[U,S,V] = svd(A - Hxeye)

%We get only the basis vectors whose singular values are sufficiently close
%to 0
mapHspace = V(:,isAlways(abs(diag(S))<svdThreshold))

%We cull the coefficients that are too small
mapHspace(abs(mapHspace) < cullingThreshold) = 0;

%which lets us get the Hamiltonian. The coefficients b are free to vary.
syms b [size(mapHspace,2) 1]
mapH = b.'*(mapHspace.'*HMonomials);
end

