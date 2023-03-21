
%% Define constants and parameters
c = 3e8; % speed of light
L = 1; % length of optical cavity
fsr = c/(2*L); % free spectral range
f = linspace(-2.5*fsr, -0.5*fsr, 1e4+1); % frequency range
r1 = 0.99; % reflectance of input mirror
r2 = 0.98; % reflectance of output mirror
t1=1-r1; t2=1-r2; 
fm = 56e6; % modulation frequency
fp=2.812e11;
w=2*pi*f;
t=L/c;
k=2*pi*f/c;
Ai=1*exp(1i*w*t);
g=1-L/r2;


%% Calculations
% Reflection coefficient of the cavity as a function of frequency
R=(sqrt(r1*r2)*(exp(1i*calculate_phi(0,0,L,k,g)-1)))./(1-r1*r2*exp(1i*calculate_phi(0,0,L,k,g)));

fn=pi*sqrt(r1*r2)/(1-r1*r2);
fprintf('The finesse is %f\n', fn);

% Calculate the Pound-Drever-Hall (PDH) readout signal by taking the product of the reflected field and the fields with the laser offset by the modulation frequency and its opposite
pdh = R.*conj(Rffm) - conj(R).*Rfnfm;

%% Plots as a function of frequency
figure();
grid on;
plot(f/1e6,abs(R), 'blue'); % Stored power
plot(f/1e6,abs(calculate_phi(0,0,L,k,g)), 'blue'); % Stored power

%% Define the phase shift in radians

function phi = calculate_phi(n, m, L, k, g)
    phi = 2*k*L-2*(n+m+1)*acos(sqrt(g));
end



