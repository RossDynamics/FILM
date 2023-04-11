function ic = ics_energy_boundary(q1,const,h,sigma,a,T,varargin)
%ICS_ENERGY_BOUNDARY Given a set of relevant variables and q1 (by default; see below), 
%calculates the rest of the initial condition (assuming 
%q_2 = p_2). If an optional variable is set to true, the function will
%instead assume p_1 was provided (in order to build out the transit picture
%when going in reverse).

if nargin == 7
    usep1 = varargin{1};
else
    usep1 = false;
end

theta = acos(a);

lambdatil = 1/T*log(sigma);

omegatil = theta / T;

p1 = const + q1;

alpha2 = 2 * h / omegatil - 2 * lambdatil / omegatil * q1 * p1;

q2 = sqrt(alpha2 / 2);

p2 = q2;

%it's kinda stupid, but if usep1 is on we just sorta assume q1 was actually
%p1 this whole time and so just interchange them
if usep1
    ic = [p1
          q2
          q1
          p2];
else    
    ic = [q1
          q2
          p1
          p2];
end

end

