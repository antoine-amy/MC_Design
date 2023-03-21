
%% Define the constants and parameters of the Fabry-Perot cavity
c = 3e8; % speed of light
L = 5e-3; % cavity length in meters
fsr = c/(2*L); % free spectral range
f = linspace(-0.5*fsr, 0.5*fsr, 1e4); % frequency range
R1 = 0.88; T1=sqrt(1-R1^2); % mirror 1
R2 = 0.99; T2=sqrt(1-R2^2); % mirror 2
k=2*pi*f/c;
g=1-L/R2;

%% Calculate the complex transmission coefficient t
ex=exp(1i*calculate_phi(0,0,L,k,g));
t=(T1*T2*ex)./(1-R1*R2*ex);

fn=pi*sqrt(R1*R2)/(1-R1*R2);
fprintf('The finesse is %f\n', fn);

%% Calculate the transmitted intensity I
I = abs(t).^2;

%% Plot the transmitted intensity as a function of wavelength
plot(f/1e6, I);
xlabel('Frequency (MHz)');
ylabel('Transmitted Intensity');
title('Transmitted Intensity of a Fabry-Perot Cavity');

%% Define the phase shift in radians

function phi = calculate_phi(n, m, L, k, g)
    phi = 2*k*L-2*(n+m+1)*acos(sqrt(g));
end

