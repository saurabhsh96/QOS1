%% SubQuestion 1; Polar and normal plots
clear;
% Giving Input Parameters
%Freq
f = 15e9;

%SoL
c = 3e8;

%Wavelength
lam = c/f;

%Dimensions of antenna
L = lam/2;
W = lam/20;
zeta = 377;

%Er
er = 1;

%Defining meshgrid
drad = pi/180;
[th, phi] = meshgrid(10e-7:drad:pi+10e-7-drad, 0:drad:2*pi-drad);
ths = 0:drad:2*pi-drad;
%Obs pt; in far field
r = 100000*lam;
[EFx, EFy, EFz] = FF(f, er, L, W, r, th, phi);
Emag = sqrt(abs(EFx).^2 + abs(EFy).^2 + abs(EFz).^2);
Emax = max(max(Emag));
%IntE = (FFTheta.^2 + FFPhi.^2);

plot(th(1,:), abs(Emag(1,:))/abs(Emax));
title('Normalized FF E-Field Intensity wrt theta, phi = 0');
ylabel('Normalized FF E-Field Intensity');
xlabel('Theta(in radians)');
ylim([0, 1]);

figure;
%polarplot(th(1,:), abs(Emag(1,:))/abs(Emax));
required = [(abs(Emag(181,:)))/(abs(Emax)) (abs(Emag(181,:)))/(abs(Emax))];
polarplot(ths, required);
title('Normalized FF E-Field Intensity wrt theta, phi = 0');

figure;
plot(th(1,:), abs(Emag(91,:))/abs(Emax))
ylim([0,1])
title('Normalized FF E-Field Intensity wrt theta, phi = 90');
ylabel('Normalized FF E-Field Intensity');
xlabel('Theta(in radians)');

figure;
required1 = [(abs(Emag(91,:)))/(abs(Emax)) (abs(Emag(271,:)))/(abs(Emax))];
polarplot(ths, required1);
title('Normalized FF E-Field Intensity wrt theta, phi = 90');

% 
% %Why normalized and real plots different
% 

%% Subquestion 2; surface plots
%Creating meshgrid for given problem, the grid takes only upper half plane
[th1, phi1] = meshgrid(10e-7:drad:pi/2+10e-7-drad, 0:drad:2*pi-drad);

Ax1 = sin(th1).*cos(phi1); %Function of sin(th).cos(phi)
Ax2 = sin(th1).*sin(phi1); %Function of sin(th).sin(phi)

EmagReq = abs(Emag(:,1:90));

figure;
surf(Ax1, Ax2, EmagReq);
title('2D plot of E-Field(magnitude) in upper medium');
zlabel('|E(th, phi)|');
xlabel('sin(th).cos(phi)');
ylabel('sin(th).sin(phi)');

%% Subquestion 3 Directivity in Broadside direction
%Define given frequency range
fRange = 5e9:1e9:30e9;

%Define dimensions of the dipole
L = 0.010; %10 mm
W = 0.001; %1 mm

%Define distance
r = 200;

%Call Directivity function %Using same meshgrid as before
DirReq = zeros(1,length(fRange));
PradReq = zeros(1,length(fRange));
a=1; %array index
for fl = fRange
    [Dir, Prad] = Directivity(fl, L, W, er, r, th, phi);
    DirReq(a) = Dir(1, 1);
    PradReq(a) = Prad;
    a=a+1;
end 

figure;
plot(fRange, DirReq);
title('Directivity');
xlabel('Frequency(in Hz)');
ylabel('D(th, phi)');