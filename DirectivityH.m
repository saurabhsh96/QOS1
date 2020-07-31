%Function to calculate the directivity of the antenna using the E-Field
function [Dir, Prad] = DirectivityH(f, L, W, er, r, thi, phi, h)
    %EfarMagSquare = Eth.^2 + Ephi.^2;
    %Finding zeta
    eps_0 = 8.854187817e-12;
    mu_0 = 1.2566370614e-6;
    zeta = round(sqrt(mu_0/(eps_0*er)));
    
    %Calculate FF
    [EFx, EFy, EFz] = FFelevated(f, er, L, W, r, thi, phi, h);
    Emag = sqrt(EFx.^2 + EFy.^2 + EFz.^2);
    
    %Calculate field intensity
    Intensity = r^2.*abs(Emag).^2./(2*zeta);
    
    %Calculate prad
    %drad is the infinitsimal change in th and phi i.e it is dth and dphi
    %drad = pi/180; 
    dth = thi(1, 2) - thi(1, 1);
    dph = phi(2, 1) - phi(1, 1);
    Prad = (sum(Intensity.*sin(thi), 'all'))*dth*dph;
    
    %Calculate directivity
    Dir = (4*pi).*Intensity/Prad; 
end