% Define parameters
c = 3e8; % speed of light
L = 0.21855230448551544; % cavity length
n = 1.1; % refractive index
R1 = 0.98; % reflectivity of mirror 1
R2 = 0.98; % reflectivity of mirror 2

% Calculate cavity finesse and free spectral range
F = pi*sqrt(R1*R2)/(1-R1*R2);
FSR = c/(2*n*L);

% Define frequency range
f = -FSR/2:FSR/100:FSR/2;

% Calculate cavity transmission for carrier wave and sidebands
Tcav_carrier = (1-R1)*(1-R2)./((1-sqrt(R1*R2))^2 + 4*sqrt(R1*R2)*sin(pi*f/FSR).^2);
Tcav_carrier = Tcav_carrier+(1-R1)*(1-R2)./((1-sqrt(R1*R2))^2 + 4*sqrt(R1*R2)*sin(pi*(f+6e6)/FSR).^2);
Tcav_carrier = Tcav_carrier+(1-R1)*(1-R2)./((1-sqrt(R1*R2))^2 + 4*sqrt(R1*R2)*sin(pi*(f+56e6)/FSR).^2);
Tcav_carrier = Tcav_carrier+(1-R1)*(1-R2)./((1-sqrt(R1*R2))^2 + 4*sqrt(R1*R2)*sin(pi*(f-6e6)/FSR).^2);
Tcav_carrier = Tcav_carrier+(1-R1)*(1-R2)./((1-sqrt(R1*R2))^2 + 4*sqrt(R1*R2)*sin(pi*(f-56e6)/FSR).^2);
% Plot transmissivity of mode cleaner cavity as a function of frequency for carrier wave and sidebands
figure;
hold on;
plot(f,Tcav_carrier,'LineWidth',3);
xlabel('Frequency (Hz)');
ylabel('Transmissivity');
title('Transmissivity of Mode Cleaner Cavity as a Function of Frequency');
legend('Carrier Wave','Sideband at 6 MHz','Sideband at 56 MHz');
hold off;