% Finding Real and Imag part of SGF with radiation condition satisfied.
clear;

%% Inputs Defination

%in Hertz (Frequency)
f = 15e9;
%in meters (speed of light)
c = 3e8; 
%in meters (Wavelegth)
lam = c/f; 
%Wavenumber (in per meter)
k0 = 2*pi/lam;

%Defining zeta
%Assuming free space
eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;
zeta = round(sqrt(mu_0/eps_0));

%Satisfying the radiation condition
obs_pt = 1000*lam;

%Defining meshgrid
drad = pi/180;
[th, phi] = meshgrid(10e-7:drad:pi+10e-7-drad, 0:drad:2*pi-drad);
% Ask, why drad and why not start from 0 also the requirement of [ths in
% instruction class?

%Defining spectral numbers
ky = 0;
kx = 0:k0/100:5*k0;

%Defining current (Elementary oriented in X)
J = [1; 0; 0]; 
%% For elementary dipole, free space components plotting

bSGF = createSGF(k0, kx, ky, zeta, 1e-7); % Is it near zero? why only one theta?
plot(kx, real(squeeze(bSGF(1,1,:,:)))); hold on; plot(kx, imag(squeeze(bSGF(1,1,:,:))));
legend('Real(SGFxx)', 'Imag(SGFxx)');
title('Real and Imag part of x-electric component of SGF wrt Kx')
xlabel('kx(rad/m)')
ylabel('Real(SGFxx) and Imag(SGFxx)')
figure;
plot(kx, real(squeeze(bSGF(3,1,:,:)))); hold on; plot(kx, imag(squeeze(bSGF(3,1,:,:))));
legend('Real(SGFzx)', 'Imag(SGFzx)');
title('Real and Imag part of z-electric field component of SGF wrt Kx')
xlabel('kz(rad/m)')
ylabel('Real(SGFzx) and Imag(SGFzx)')