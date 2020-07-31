% This code is for writing genaralized SGF
function SGF = createSGF(k0, kx, ky, zeta, th)
    kz = -1j*sqrt(-(k0.^2-kx.^2-ky.^2));
    %kz = sqrt(k0^2-kx.^2-ky.^2);
    C = -zeta./(2.*k0.*kz);
    if length(th)==1 
        if th <= pi/2
            Dxx = k0.^2 - kx.^2;  Dxy = -kx.*ky;         Dxz = -kx.*kz;
            Dyx = -ky.*kx;       Dyy = k0^2 - ky.^2;     Dyz = -ky.*kz;
            Dzx = -kz.*kx;       Dzy = -kz.*ky;          Dzz = k0^2 - kz.^2;
        else
            Dxx = k0^2 - kx.^2;  Dxy = -kx.*ky;          Dxz = kx.*kz;
            Dyx = -ky.*kx;       Dyy = k0^2 - ky.^2;     Dyz = ky.*kz;
            Dzx = kz.*kx;        Dzy = kz.*ky;           Dzz = k0^2 - kz.^2;
        end
    else
        %Defining Length of theta
        thetaL = length(th(1,:));
        Dxx = k0.^2 - kx.^2;  
        Dxy = -kx.*ky;          
        Dxz = -kx.*kz;       Dxz(:,1+thetaL/2:thetaL) = -Dxz(:,1+thetaL/2:thetaL);
        Dyx = -ky.*kx;       
        Dyy = k0.^2 - ky.^2;     
        Dyz = -ky.*kz;       Dyz(:,1+thetaL/2:thetaL) = -Dyz(:,1+thetaL/2:thetaL);
        Dzx = -kz.*kx;       Dzx(:,1+thetaL/2:thetaL) = -Dzx(:,1+thetaL/2:thetaL);
        Dzy = -kz.*ky;       Dzy(:,1+thetaL/2:thetaL) = -Dzy(:,1+thetaL/2:thetaL);
        Dzz = k0.^2 - kz.^2;
    end
    SGF(1,1,:,:) = C.*Dxx; SGF(1,2,:,:) = C.*Dxy; SGF(1,3,:,:) = C.*Dxz;
    SGF(2,1,:,:) = C.*Dyx; SGF(2,2,:,:) = C.*Dyy; SGF(2,3,:,:) = C.*Dyz;
    SGF(3,1,:,:) = C.*Dzx; SGF(3,2,:,:) = C.*Dzy; SGF(3,3,:,:) = C.*Dzz;
end
