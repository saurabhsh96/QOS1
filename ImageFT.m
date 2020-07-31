%Modified FT of current for PEC plane
function [currFT] = ImageFT(k0, kl, kw, L, W, J, kh, h)
    %kl represents the orientation of Length of dipole
    %kw represents the orientation of Width of dipole
    %Assuming PWS as current Dist.
    %Assuming I0 = 1
    %J = [1; 0; 0];
    %h = distance above the plane
    %kh = k along the height axis of the dipole
    Tx = sinc(kw*W/2/pi); %Why divide by pi?
    Lx = (4*1j*k0*(cos(kl*L/2) - cos(k0*L/2)).*sin(kh*h)./((k0.^2 - kl.^2).*sin(k0*L/2)));
    %FT of delta(z-h)-delta(z+h) is 2jsin(kz*h)
    currFT(1,:,:) = Lx.*Tx*J(1);
    currFT(2,:,:) = Lx.*Tx*J(2);
    currFT(3,:,:) = Lx.*Tx*J(3);
end