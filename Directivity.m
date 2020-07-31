%Function to calculate the directivity of the antenna using the E-Field
function [Dir, Prad] = Directivity(f, L, W, er, r, thi, phi)
    %EfarMagSquare = Eth.^2 + Ephi.^2;
    %Finding zeta
%     eps_0 = 8.854187817e-12;
%     mu_0 = 1.2566370614e-6;
%     zeta = (sqrt(mu_0/(eps_0*er)));
    zeta0 = 120*pi;
    zetaS = zeta0./sqrt(er);
    
    %Calculate FF
    [EFx, EFy, EFz] = FF(f, er, L, W, r, thi, phi);
    Emag = sqrt((abs(EFx)).^2 + (abs(EFy)).^2 + (abs(EFz)).^2);
    
    %Calculate field intensity
    Intensity = r^2.*((abs(Emag)).^2)./(2*zetaS);
    
    %Calculate prad
    %drad is the infinitsimal change in th and phi i.e it is dth and dphi
    %drad = pi/180; 
    dth = thi(1, 2) - thi(1, 1);
    dph = phi(2, 1) - phi(1, 1);
    Prad = (sum(Intensity.*sin(thi), 'all'))*dth*dph;
    
    %Calculate directivity
    Dir = (4*pi).*Intensity/Prad; 
end