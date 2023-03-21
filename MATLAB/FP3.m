% Define parameters
c = 3e8; % speed of light
F=200; %finesse of the cavity
lbda=1064;

L = 3.6e-2; % cavity length in meters
fsr = c/(2*L); % free spectral range
f = linspace(-0.5*fsr, 0.5*fsr, 1e4); % frequency range

r1 = 0.939; T1=sqrt(1-r1^2); % mirror 1
r2 = 0.939; T2=sqrt(1-r2^2); % mirror 2
n_modes = 7; % Total number of modes (n+m)
k=2*pi*f/c;


%% Calculate the complex transmission coefficient t
ex=exp(1i*calculate_phi(0,0,L,k,g));
t=(T1*T2*ex)./(1-r1*r2*ex);

% Calculate round-trip phase shift
phi = 4*pi*L;

% Create input field with n modes
E_in = ones(n_modes,1);

% Calculate transmission coefficients for each mode
T = zeros(n_modes);
for i = 1:n_modes
    %T(i) = sqrt(r1)*sqrt(r2)*exp(1i*(ii-jj)*phi)/(1-r1*r2*exp(1i*(ii-jj)*phi));
    T(i) = 1/(1+(2*F/pi)^2*sin(i*acos(sqrt(g))+pi/2)^2);
end

% Calculate output field with m modes
E_out = T.*E_in;

% Plot output field modes intensity
figure;
stem(abs(E_out));
title(sprintf('Output field with %d modes intensity',n_modes-2));

%% Define the phase shift in radians

function phi = calculate_phi(n, m, L, k, g)
    phi = 2*k*L-2*(n+m+1)*acos(sqrt(g));
end