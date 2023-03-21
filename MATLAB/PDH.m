%% Define constants and parameters
c = 3e8; % speed of light
L = 1; % length of optical cavity
fsr = c/(2*L); % free spectral range
f = linspace(-0.5*fsr, 0.5*fsr, 1e4+1); % frequency range
f = f(1:end-1); % remove last element to match size of other variables
r1 = 0.99; % reflectance of input mirror
r2 = 0.98; % reflectance of output mirror
t1 = sqrt(1-r1^2); % transmissivity of input mirror
fm = 56e6; % modulation frequency



%% Calculations
% Reflection coefficient of the cavity as a function of frequency
R=(r1-(r1^2+t1^2)*r2*exp(2i*pi*f*L/c))./(1-r1*r2*exp(2i*pi*f*L/c));
Rffm=(r1 - (r1^2 + t1^2)*r2*exp(2i*pi*(f + fm)*L/c)) ./ (1 - r1*r2*exp(2i*pi*(f + fm)*L/c)); % when the laser is offset by the modulation frequency
Rfnfm = (r1 - (r1^2 + t1^2)*r2*exp(2i*pi*(f - fm)*L/c)) ./ (1 - r1*r2*exp(2i*pi*(f - fm)*L/c)); % when the laser is offset by the opposite of the modulation frequency

fn=pi * sqrt(r1*r2)/(1-r1*r2);
fprintf('The finesse is %f\n', fn);

% Calculate the Pound-Drever-Hall (PDH) readout signal by taking the product of the reflected field and the fields with the laser offset by the modulation frequency and its opposite
pdh = R.*conj(Rffm) - conj(R).*Rfnfm;

%% Plots as a function of frequency
figure();
grid on;
yyaxis left;
plot(f/1e6, 100*abs(R).^2, 'blue'); % Stored power
ylabel('Reflected power (%)');
ylim([-100, 100]);
hold on;
yyaxis right;
plot(f/1e6, 180*angle(R)/pi, 'red');  % Reflected phase
ylabel('Reflected phase (deg.)');
ylim([-100, 100]);
xlabel('Frequency (MHz)');
plot(f/1e6, 100*imag(pdh), 'black'); % PDH readout
title('PDH technique in a FP cavity');
legend('Stored power', 'Reflected phase', 'PDH readout (arb.)', 'Location', 'northwest');



