function ic = ics_energy_boundaryER3BP(q1,const,h,sigma,a,T)
%ICS_ALONG_LINE Given a set of relevant variables and q1, 
%calculates the rest of the initial condition (assuming 
%q_2 = p_2). Adds on an initial condition for the true anomaly.

theta = acos(a);

lambdatil = 1/T*log(sigma);

omegatil = theta / T;

p1 = const + q1;

alpha2 = 2 * h / omegatil - 2 * lambdatil / omegatil * q1 * p1;

q2 = sqrt(alpha2 / 2);

p2 = q2;

ic = [q1
      q2
      p1
      p2
      0];

end

