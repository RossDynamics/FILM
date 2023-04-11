function opt = multiRefineHandle(p1val,q2,p2,h,ptMaker,H2energy,revquadrat)
    %Used to refine the correct q2 and p2 to make the energy and center
    %angle match the correct values after an iterate.
    pt = revquadrat(ptMaker(p1val,q2,p2));
    energyopt = real(log(H2energy(0,pt)) - log(h)); %smootherabs(H2energy(0,pt) - h,1e+15)/h;
    angleopt = (pt(2) - pt(4))/abs(pt(2)); %smootherabs(pt(2) - pt(4),1e+15);
      
    disp('criteria:')
    disp(energyopt)
    disp(angleopt)
    disp('H2:')
    disp(H2energy(0,pt))
    disp('pt(2):')
    disp(pt(2))
    
    opt = [energyopt
           angleopt];
end

