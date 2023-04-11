%Uses a "smarter" version of the approach in analyzeTransitQuadratic to
%refine the transit/nontransit boundary into a quadratic approximation of
%the stable (or unstable) manifold. Run analyzeMonodromy or 
%analyzeMonodromyER3BP first.

%%
%PARAMETERS
%Set equal to 1 for the stable manifold and -1 for the unstable manifold
direction = -1;

%Change depending on the model:
%period = getSolarPeriod(c);
period = 2*pi;

%The array of p1's (for the stable manifold) or q1's (for the stable
%manifold) at which the manifold will be sampled.
samplePts = linspace(mp('-1e-5'), mp('1e-5'), 100);

%The distance between transit orbits and nontransit orbits must be below
%this tolerance for the solver to stop at a given sample point.
distTol = mp('1e-16');

%The quadratic energy
h = mp('1e-14');

%The range within which we are searching for the manifold; \pm initialGuess
%are used as the initial bounds for the bisection-ish algorithm we are 
%employing.
initialGuess = mp('1e-8');

%We need these for the energy refiner. The underlying values should have
%already been computed by analyzeMonodromy for you, so you shouldn't need
%to mess with these.
tillambda = log(cg(c,'p.sigma'))/period;
tilomega = acos(cg(c,'p.a'))/period;

%%
%THE QUADRATIC MAP

%The initial guess for the transit orbit will 

c = coordreset(c);
c = useMomentum(c);

%We get the quadratic map. It is faster if we try to compute it within the
%standard momentum frame and then include conversion functionality manually
% disp('Calculating state transition tensors:')
% [phi,~,~] = stm(linspace(0,direction * period,3),mp(cg(c,'lm.y0')),c,[],2);
% disp('Calculation complete.')

eigp = eigenbasis.basis.value;
eigo = eigenbasis.origin.value;

quadratstd = @(x)(phi{1}(1:4,1:4,end)*x +...
    getfield(1/2*ttv(ttv(tensor(phi{2}(1:4,1:4,1:4,end)),x,2),x,2),'data'));
quadrat=@(x) eigp \ quadratstd(eigp * x);

linstd = @(x)phi{1}(1:4,1:4,end)*x;
lin=@(x) eigp \ linstd(eigp * x);

%Let's set the eigenbasis now
c = coordset(c,eigenbasis);

%%
%ENERGY FIXER

%For some q1 and p1, we calculate the q2 = p2 necessary in the center
%projection such that the initial condition has quadratic energy h.
refineq2 = @(q1,p1) sqrt((h - tillambda*q1*p1)/tilomega);

%From there, we can build out a function handle to make the ic. I know
%there's an unnecessary function call to refineRad I know it could be
%more optimized don't @ me
getIC = @(q1,p1)[q1 refineq2(q1,p1) p1 -refineq2(q1,p1)].'

%We create a wrapper function so that, instead of inputing q1 and p1,
%directly, the bisection solver sees slots for the *sampled and refined*
%variables. The sampled variable is the fixed one drawn from samplePts; the
%refined variable is the one to be refined. That is, we leave or
%interchange q1 and p1 depending upon their roles in the solver so that the
%refiner doesn't have to deal with that specific detail.
switch direction
    case 1
        getICrs = @(samp,refnd)getIC(refnd,samp);
    case -1
        getICrs = @(samp,refnd)getIC(samp,refnd);
end

%%
%REFINEMENT

manifold = NaN(4,numel(samplePts));
for i = 1:numel(samplePts)
     lowerBnd = -initialGuess;
     upperBnd = initialGuess;
     
     %We refine by finding the smallest range of values such that one
     %transits and one does not. We continue this process until the
     %distance between the bounds is sufficiently small.
     while abs(upperBnd - lowerBnd) >= distTol
         
         %We determine whether the midpoint of the current bounds is
         %transit or nontransit
         midVal = (upperBnd + lowerBnd)/2;
         testPt = getICrs(samplePts(i),midVal);
         
%          if ~isreal(testPt)
%              upperBnd = NaN;
%              disp('Complex test point detected. Skipping iterate.')
%              continue
%          end
%          
         %We iterate the test point forward one iteration under the
         %quadratic map (which could actually be backwards, depending on
         %which manifold you're targeting)
         testPtFd = quadrat(testPt);
         
%          c = cs(c,'s.o.v.dmode','1');
%          cplot(testPtFd,c,'ok');
%          drawnow;
%          c = cs(c,'s.o.v.dmode','position');
         
         %Transit evidence is the value of the iterated coordinate relevant
         %to distinguishing transit from nontransit behavior. For example, 
         %if we are refining the stable manifold, it is the iterated value
         %of q1.
         switch direction
             case 1
                 transitEvidence = testPtFd(1);
             case -1
                 transitEvidence = testPtFd(3);
         end
         
         %Now, we need to figure out whether the transit evidence means
         %that the point is transit or nontransit. This occurs
         %under the following condition:
         %isTransit = transitEvidence * sign(samplePts) > 0;
         %Then, as the lower bound is nontransit in the upper halfplane
         %and transit in the lower halfplane, whereas the upper bound is
         %transit in the former and nontransit in the lower, we can choose
         %which bound midVal should become. Specifically, it should replace
         %the existing bound that it matches.
         %However, there is a simpler approach. The upper bound should
         %always go to infinity under repeated iterates and the lower bound
         %should always go to negative infinity. Thus, the following also
         %works:
         if transitEvidence > 0
             upperBnd = midVal;
         elseif transitEvidence < 0
             lowerBnd = midVal;
         else
             %If this happens, we have somehow found a point perfectly on
             %the stable manifold (I know not how), so we just break.
             break;
         end
         
     end
     
     manifold(:,i) = getICrs(samplePts(i),(upperBnd + lowerBnd)/2);
     disp(i)
     
end

 c = cs(c,'s.o.v.dmode','1');
 cplot(manifold,c,'-k');
 drawnow;
 c = cs(c,'s.o.v.dmode','position');