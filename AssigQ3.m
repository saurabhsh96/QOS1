%% SubQuestion 1; Polar and normal plots
clear;
% Giving Input Parameters
%Freq
%f = 15e9;
f = 10e9;
%SoL
c = 3e8;

%Wavelength
lam = c/f;

%Dimensions of antenna
% L = lam/2;
% W = lam/20;
L = lam/200000;
W = lam/200000;

zeta = 377;

%Elevation
%h = 0.015;
h = 0.002;

%Er
er = 1;

%Defining meshgrid
drad = pi/180;
[th, phi] = meshgrid(10e-7:drad:pi+10e-7-drad, 0:drad:2*pi-drad);
ths = 0:drad:2*pi-drad;

%Obs pt; in far field
%r = 100000*lam;
r = h;
% Image Current
[EFx, EFy, EFz] = FFelevated(f, er, L, W, r, th, phi, h);
Emag = sqrt(EFx.^2 + EFy.^2 + EFz.^2);
%Emag(:,91:180) = 0;
Emax = max(max(Emag));

% Without PEC
[EFx1, EFy1, EFz1] = FF(f, er, L, W, r, th, phi);
EmagWP = sqrt(EFx1.^2 + EFy1.^2 + EFz1.^2);
EmaxWP = max(max(EmagWP)); %Technically multiplocation by e^jkzh has to be done to J

%IntE = (FFTheta.^2 + FFPhi.^2);
%Why not getting correct for 1:180 plot itself?
Emag1 = Emag;
Emag1(:,91:180) = 0;
plot(th(1,:), abs(Emag1(1,:))/abs(EmaxWP)); hold on;
plot(th(1,:), abs(EmagWP(1,:))/abs(EmaxWP));
title('Normalized FF E-Field Intensity wrt theta, phi = 0');
ylabel('Normalized FF E-Field Intensity');
xlabel('Theta(in radians)');
legend('With PEC', 'Without PEC');
ylim([0, 2]);

figure;
Emag2=Emag;
Emag2(:,1:90) = 0;
required = [(abs(Emag1(1,:)))/(abs(EmaxWP)) (abs(Emag2(1,:)))/(abs(EmaxWP))];
requiredWP = [(abs(EmagWP(181,:)))/(abs(EmaxWP)) (abs(EmagWP(181,:)))/(abs(EmaxWP))];
%required = [(abs(Emag(1,:)))/(abs(Emax)) (abs(Emag(1,:)))/(abs(Emax))];
%polarplot(th(1,:), abs(Emag2(1,:))/(abs(Emax)));
polarplot(ths, required); hold on;
polarplot(ths, requiredWP);
title('Normalized FF E-Field Intensity wrt theta, phi = 0');
legend('With PEC', 'Without PEC');

figure;
plot(th(1,:), abs(Emag1(91,:))/abs(EmaxWP)); hold on;
plot(th(1,:), abs(EmagWP(91,:))/abs(EmaxWP));
ylim([0,2])
title('Normalized FF E-Field Intensity wrt theta, phi = 90');
ylabel('Normalized FF E-Field Intensity');
xlabel('Theta(in radians)');
legend('With PEC', 'Without PEC');

figure;
required1WP = [(abs(EmagWP(91,:)))/(abs(EmaxWP)) (abs(EmagWP(271,:)))/(abs(EmaxWP))];
required1 = [(abs(Emag1(91,:)))/(abs(EmaxWP)) (abs(Emag2(91,:)))/(abs(EmaxWP))];
polarplot(ths, required1); hold on;
polarplot(ths, required1WP);
title('Normalized FF E-Field Intensity wrt theta, phi = 90');
legend('With PEC', 'Without PEC');

%% Subquestion 2 Directivity

%Define given height range
hRange = 1e-3:1e-4:20e-3;

%Call Directivity function %Using same meshgrid as before
DirReq = zeros(1,length(hRange));
PradReq = zeros(1,length(hRange));
a=1; %array index
for hl = hRange
    [Dir, Prad] = DirectivityH(f, L, W, er, r, th, phi, hl);
    DirReq(a) = Dir(1, 1);
    PradReq(a) = Prad;
    a=a+1;
end 

figure;
plot(hRange.*10^3, DirReq);
title('Directivity');
xlabel('Height(in mm)');
ylabel('D(th, phi)');