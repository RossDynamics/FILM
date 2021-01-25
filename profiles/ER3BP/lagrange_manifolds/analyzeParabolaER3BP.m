%Used for finding the bounding parabola using the symplectic eigenbasis.
%Run analyzeMonodromyER3BP first.

%Let's set the eigenbasis first
c = coordset(c,eigenbasis);

parabola = [];

c = cs(c,'s.o.v.dmode','1');

for i = logspace(-9,-2,5000)

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

    %Now, we build the array using ics_energy_boundary:
    arrayics = ics_energy_boundaryER3BP(q1,const,h,...
                                   cg(c,'p.sigma'),cg(c,'p.a'),...
                                   cg(c,'p.T'),newy0(5));

    while true
        q1 = q1 + arrayStep;
        ic = ics_energy_boundaryER3BP(q1,const,h,...
                                   cg(c,'p.sigma'),cg(c,'p.a'),...
                                   cg(c,'p.T'),newy0(5));
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