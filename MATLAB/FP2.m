
%% Define the constants and parameters of the Fabry-Perot cavity
c = 3e8; % speed of light
lambda = input("Enter the desired wavelength (in nm): ")*1e-9;
R1 = 0.939; T1=sqrt(1-R1^2); % mirror 1
R2 = 0.939; T2=sqrt(1-R2^2); % mirror 2
g=0.9; q=6532;

% Calculate the necessary length of the cavity

L =(lambda/2)*(q+0.5+(acos(sqrt(g))/pi));
fsr = c/(2*L); % free spectral range
f = linspace(-0.5*fsr, 0.5*fsr, 1e4); % frequency range
k=2*pi*f/c;

%% Calculate the complex transmission coefficient t
ex=exp(1i*calculate_phi(0,0,L,k,g));
t=(T1*T2*ex)./(1-R1*R2*ex);

fn=pi*sqrt(R1*R2)/(1-R1*R2);
fprintf('Finesse: %f\n', fn);
Pmax=(T1/(1-R1*R2))^2;
fprintf('Cavity lenght (mm): %f\n', L*1e3);

%% Calculate the transmitted intensity I
I = abs(t).^2;

%% Plot the transmitted intensity as a function of wavelength
plot(f/1e6, I);
xlabel('Frequency (MHz)');
ylabel('Transmitted Intensity');
title(sprintf('Transmitted Intensity of a Fabry-Perot Cavity at %g nm', lambda*1e9));

%% Define the phase shift in radians

function phi = calculate_phi(n, m, L, k, g)
    phi = 2*k*L-2*(n+m+1)*acos(sqrt(g));
end

