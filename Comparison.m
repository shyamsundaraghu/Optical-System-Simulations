%% Chirped Fiber Bragg Grating (CFBG) Simulation
clc; clear; close all;

%% Parameters
L = 10e-3;               % Grating length [m]
neff = 1.45;             % Effective refractive index
Delta_n = 1e-4;          % Refractive index modulation
Lambda0 = 1550e-9;       % Start period [m]
c = 3e8;                 % Speed of light [m/s]
x = 1e-6;                % Linear chirp coefficient [m/m]
numPoints = 1000;        % Number of points along grating

%% Position along grating
z = linspace(0, L, numPoints);

%% Grating period (linear chirp)
Lambda = Lambda0 + x*z;

%% Apodization functions
% Blackman
Xb = 2*pi*(z - L/2)/L;
g_blackman = (1 + 1.19*cos(Xb) + 0.19*cos(2*Xb))/2.38;

% Cauchy
C = 0.5;
g_cauchy = 1 - ((z - L/2)/(L/2)).^2 - 2*C*((z - L/2)/(L/2)).^2;

% Hamming
H = 0.9;
g_hamming = (1 + H*cos(2*pi*(z - L/2)/L))/(1 + H);

% Gaussian
a = 1;
g_gaussian = exp(-a*((z/L - 0.5).^2));

% Sinc
g_sinc = sin(6*(z/L - 0.5))./(6*(z/L - 0.5));

% Hyper-tangent
s = 6;
g_tanh = 0.5*(tanh(s*(z/L - 0.5)) + 1);

%% Effective refractive index along grating
n_z_blackman = neff + Delta_n*g_blackman.*cos(2*pi*z./Lambda);
n_z_cauchy    = neff + Delta_n*g_cauchy.*cos(2*pi*z./Lambda);
n_z_hamming   = neff + Delta_n*g_hamming.*cos(2*pi*z./Lambda);
n_z_gaussian  = neff + Delta_n*g_gaussian.*cos(2*pi*z./Lambda);
n_z_sinc      = neff + Delta_n*g_sinc.*cos(2*pi*z./Lambda);
n_z_tanh      = neff + Delta_n*g_tanh.*cos(2*pi*z./Lambda);

%% Bragg wavelength along the grating
lambdaB_blackman = 2*n_z_blackman.*Lambda;
lambdaB_cauchy    = 2*n_z_cauchy.*Lambda;
lambdaB_hamming   = 2*n_z_hamming.*Lambda;
lambdaB_gaussian  = 2*n_z_gaussian.*Lambda;
lambdaB_sinc      = 2*n_z_sinc.*Lambda;
lambdaB_tanh      = 2*n_z_tanh.*Lambda;

%% Chirp range
Delta_lambda = 2*neff*(max(Lambda) - min(Lambda));

%% Time delay for each wavelength
tau_blackman = ((lambdaB_blackman - Lambda0)./(2*neff*Delta_lambda)) * L/c;

%% Plotting
figure;
plot(z*1e3, lambdaB_blackman*1e9, 'r', 'LineWidth', 1.5); hold on;
plot(z*1e3, lambdaB_gaussian*1e9, 'b', 'LineWidth', 1.5);
plot(z*1e3, lambdaB_hamming*1e9, 'g', 'LineWidth', 1.5);
xlabel('Position along grating z [mm]');
ylabel('Bragg wavelength \lambda_B [nm]');
legend('Blackman', 'Gaussian', 'Hamming');
title('Bragg wavelength along Chirped Fiber Bragg Grating');
grid on;

figure;
plot(z*1e3, g_blackman, 'r', 'LineWidth', 1.5); hold on;
plot(z*1e3, g_gaussian, 'b', 'LineWidth', 1.5);
plot(z*1e3, g_hamming, 'g', 'LineWidth', 1.5);
xlabel('Position along grating z [mm]');
ylabel('Apodization function g(z)');
legend('Blackman', 'Gaussian', 'Hamming');
title('Apodization functions along the grating');
grid on;

disp(['Chirp range Delta lambda: ', num2str(Delta_lambda*1e9), ' nm']);
