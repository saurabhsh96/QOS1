% Routine to calculate FF using SGF and FT
function [EFx, EFy, EFz] = FF(f, er, L, W, r, thi, phi)
    %Calculating required parameters from given
    %Speed of light
    c = 3e8; 
    
    %Wavelength (in m)
    lam = c/f;
    k0 = 2*pi/lam;
    
    %Impedance (in Ohm)
    eps_0 = 8.854187817e-12;
    mu_0 = 1.2566370614e-6;
    zeta = round(sqrt(mu_0/(eps_0*er)));
    
    %Defining spectral numbers
    kx = k0.*sin(thi).*cos(phi);
    ky = k0.*sin(thi).*sin(phi);
    kz = k0.*cos(thi);
    %kz = sqrt(k0^2-kx.^2-ky.^2);
    
    %Defining current distribution %This has to be changed to get a better
    %and more inclusive function, it should be an input, anyway, for now
    %keep it as it is
    J = [1, 0, 0];
    
    %Calling FT J oriented along Kx in length, so, kl = kx and kw = ky
    %If the source is elementary then this changes
    if (L==0 && W==0)
        FFCurrentFT = ones([size(J, 2) size(kx)]);
        FFCurrentFT(2, :, :) = 0;
        FFCurrentFT(3, :, :) = 0;
    else 
        FFCurrentFT = CurrentFT(k0, kx, ky, L, W, J);
    end 
    
    %Calling SGF
    FFSGF = createSGF(k0, kx, ky, zeta, 1);
    
    %Finding Electric field; Assuming obs pi in far field
    EFx = 1j.*kz.*(squeeze(FFSGF(1,1,:,:)).*squeeze(FFCurrentFT(1,:,:)) ...
        + squeeze(FFSGF(1,2,:,:)).*squeeze(FFCurrentFT(2,:,:)) ...
        + squeeze(FFSGF(1,3,:,:)).*squeeze(FFCurrentFT(3,:,:))).*(exp(-1j*k0*r)/(2*pi*r));
    
    EFy = 1j.*kz.*(squeeze(FFSGF(2,1,:,:)).*squeeze(FFCurrentFT(1,:,:)) ...
        + squeeze(FFSGF(2,2,:,:)).*squeeze(FFCurrentFT(2,:,:)) ...
        + squeeze(FFSGF(2,3,:,:)).*squeeze(FFCurrentFT(3,:,:))).*(exp(-1j*k0*r)/(2*pi*r));
    
    EFz = 1j.*kz.*(squeeze(FFSGF(3,1,:,:)).*squeeze(FFCurrentFT(1,:,:)) ...
        + squeeze(FFSGF(3,2,:,:)).*squeeze(FFCurrentFT(2,:,:)) ...
        + squeeze(FFSGF(3,3,:,:)).*squeeze(FFCurrentFT(3,:,:))).*(exp(-1j*k0*r)/(2*pi*r));
    
    
    %FFR = EFx.*sin(th).*cos(phi) + EFy.*sin(th).*sin(phi) - EFz.*cos(th);
    %FFTheta = EFx.*cos(th).*cos(phi) + EFy.*cos(th).*sin(phi) - EFz.*sin(th);
    %FFTheta = EFy.*cos(th) - EFz.*sin(th); 
    %FFPhi = -EFx.*sin(phi) + EFy.*cos(phi);
    
    %Assuming obs pt is in far field
%     EF = zeros(3,360,180);
%     for k=1:360
%         for m=1:180
%             EF(:,k,m) = 1j*kz(k,m)*FFSGF(:,:,k,m)*FFCurrentFT(:,k,m)*(exp(-1j*k0*r)/(2*pi*r));
%         end
%     end
%  
%       %Theta and Phi components of Electric Field
%       FFTheta = zeros(360,180);
%       FFPhi = zeros(360,180);
%       for k=1:360
%         for m=1:180
%             %EF(:,k,m) = 1j*kz(k,m)*FFSGF(:,:,k,m)*FFCurrentFT(:,k,m)*(exp(-1j*k0*r)/(2*pi*r));
%             FFTheta(k,m) = EF(1,k,m).*cos(th(k,m)).*cos(phi(k,m)) + EF(2,k,m).*cos(th(k,m)).*sin(phi(k,m)) - EF(3,k,m).*cos(th(k,m));
%             FFPhi(k,m) = -EF(1,k,m).*sin(phi(k,m)) + EF(2,k,m).*cos(phi(k,m));
%         end
%       end 
end
