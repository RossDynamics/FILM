%Used for finding the bounding parabola using the symplectic eigenbasis.
%Run analyzeMonodromy first.

%Let's set the eigenbasis first
c = coordset(c,eigenbasis);

%To consider initial conditions in the upper halfplane, set halfplane equal
%to 1. To consider initial conditions in the lower halfplane, set
%halfplane equal to -1.
halfplane = 1;

%We consider small h = H2(x).
h = 1e-7;

parabola = [];

c = cs(c,'s.o.v.dmode','1');

for i = logspace(-9,-2,200)

    %We create an array of initial conditions from within the
    %symplectic eigenbasis along the bounding line. We set const (where p1 =
    %const + q1) small

    const = halfplane * i;

    %We also set a beginning for the array's q1's. This beginning coincides 
    %with the p1 + q1 = 0 line, so we have q1 = -const / 2.

    q1 = -const / 2;

    %The other end of
    %the array will be found dynamically; it will occur when the initial
    %conditions become complex. We specify a small arrayStep for creating 
    %subsequent q2's.

    arrayStep = halfplane * 5e-7;

    %Now, we build the array using ics_energy_boundary:
    arrayics = ics_energy_boundary(q1,const,h,...
                                   cg(c,'p.sigma'),cg(c,'p.a'),cg(c,'p.T'));

    while true
        q1 = q1 + arrayStep;
        ic = ics_energy_boundary(q1,const,h,...
                                   cg(c,'p.sigma'),cg(c,'p.a'),cg(c,'p.T'));
        %If we get complex values, we know we've reached the boundary
        if ~isreal(ic)
            break;
        end
        arrayics = [arrayics ic];
    end
    
    %cplot(arrayics,c,'ok');
    
    parabola = [parabola arrayics(:,end)];
end

cplot(parabola,c,'-k');

c = cs(c,'s.o.v.dmode','position');