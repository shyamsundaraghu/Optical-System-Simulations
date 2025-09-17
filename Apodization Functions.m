% cfbg_sim.m
% Chirped Fiber Bragg Grating simulation using coupled-mode theory (RK4)
% Supports several apodization profiles (Blackman, Cauchy, Hamming, Gaussian, Sinc, tanh, Sinc Cauchy)
% Outputs reflectivity vs wavelength and shows apodization along grating.

clear; ; clc;

%%%%%% User parameters (change as needed) %%%%%%
L = 50e-3;                % grating length (m) -> 50 mm from paper
neff = 1.45;              % effective index
delta_n_peak = 1e-4;      % peak index modulation (Δn)
lambda0 = 1550e-9;        % central wavelength (m)
chirp_delta_lambda = 0.25e-9; % total chirp across grating in wavelength (m) (paper: 0.25 nm)
Nz = 4000;                % number of z-slices (increase for accuracy)
wvl_span = 1546e-9:0.01e-9:1554e-9;  % wavelength vector to sweep (m)
apod_type = 'ellipse';   % 'blackman','cauchy','hamming','gaussian','sinc','tanh','proposed'
% Apodization params (paper-inspired defaults)
A_sinc = 1; C_cauchy = 0.3; s_param = 6; H = 0.9; a_gauss = 1;

%%%% Derived parameters %%%%
z = linspace(0,L,Nz);
dz = z(2)-z(1);

% Choose linear chirp in grating period Lambda(z)
% Base period Lambda0 set so that central Bragg wavelength (lambda0) corresponds to Lambda0
Lambda0 = lambda0/(2*neff);
% Convert chirp in wavelength to chirp in Lambda:
% Δλ_total = 2*neff*(Λ_long - Λ_short)  => ΔΛ = Δλ_total/(2*neff)
DeltaLambda_total = chirp_delta_lambda/(2*neff);
Lambda_z = Lambda0 + ( ( (0:Nz-1)/(Nz-1) ) - 0.0 ) * DeltaLambda_total;  % linear chirp across length
Lambda_z = Lambda0 + (z/L)*DeltaLambda_total;  % nicer mapping

% Local Bragg wavelength along grating
lambdaB_z = 2*neff.*Lambda_z;

% apodization profiles (normalized between 0 and 1)
u = (z - L/2)./L;           % normalized coordinate: -0.5..+0.5 -> u in [-0.5,0.5]
x = 2*pi*u;                % generic variable used in paper simplifications

switch apod_type
    case 'blackman'
        % Blackman as in paper eq (7)-(8) [normalized]
        g = (1 + 1.19*cos(x) + 0.19*cos(2*x))/2.38;
        g(g<0) = 0;
    case 'cauchy'
        % simplified Cauchy-like used in paper (Eq.9 simplified)
        % use C parameter
        C = 0.5;
        g = 1 - ( (2*u).^2 )./( 1 + 2*C*(2*u).^2 );  % a smooth shape, positive
        g = (g - min(g))/(max(g)-min(g)); % normalize
    case 'hamming'
        g = (1 + H*cos(2*pi*(z - L/2)/L))./(1 + H);
        g = (g - min(g))/(max(g)-min(g));
    case 'gaussian'
        g = exp(-a_gauss * ( (z./L) - 0.5 ).^2 );
        g = (g - min(g))/(max(g)-min(g));
    case 'sinc'
        y = 2*pi*(z - L/2)/L;     % range roughly -pi..pi
        y(y==0) = eps;
        A = 3;
        g = A * sin(y)./y;
        g = (g - min(g))/(max(g)-min(g));
    case 'tanh'
        s = s_param;
        g = tanh(s*z/L)./tanh(s) .* tanh(s*(1 - z./L))./tanh(s);
        g = (g - min(g))/(max(g)-min(g));
    case 'Uniform'
        g=1;
    case 'Sinc Cauchy'

        xnorm = 2*(z - L/2)/L;
        % simplified Cauchy-like:
        C = C_cauchy; % 0.3 recommended in paper
        g_c = (1 - xnorm.^2)./(1 - (C*xnorm).^2 + 1e-12);
        % sinc part
        A = A_sinc;
        y = pi * xnorm;  % maps -1..1 -> -pi..pi
        y(y==0) = eps;
        s_f = A .* sin(y)./y;
        % product
        g = g_c .* s_f;
        % normalize to [0,1]
        g = (g - min(g))/(max(g)-min(g));
    otherwise
        error('Unknown apodization type.');
end

% coupling coefficient kappa(z). Approx formula: kappa = pi * delta_n(z) / lambda
% delta_n(z) = delta_n_peak * g(z)
kappa_z = pi * delta_n_peak .* g ./ lambda0;   % (1/m)
% For more accurate modeling you can include mode overlap factor. Here assumed =1

% Preallocate reflectivity array
R = zeros(size(wvl_span));

% Numerically integrate coupled-mode equations using RK4 for each lambda
% Forward A(z) and backward B(z) amplitudes. Standard CMT:
% dA/dz =  j * kappa(z) * B .* exp(+2j * delta(z) .* z)
% dB/dz = -j * kappa(z) * A .* exp(-2j * delta(z) .* z)
% where delta(z) = (beta - pi/Lambda(z)) / 2, beta = 2*pi*neff/lambda
% We'll evaluate delta(z) per wavelength and per z.

fprintf('Simulating %d wavelengths, Nz = %d ...\n', length(wvl_span), Nz);

for iw = 1:length(wvl_span)
    lambda = wvl_span(iw);
    beta = 2*pi*neff / lambda;                 % propagation constant
    % local Bragg propagation constant = pi / Lambda(z)  (since Bragg condition beta=pi/Lambda)
    delta_z = (beta - pi./Lambda_z)/2;         % detuning function [1/m]
    
    % initial conditions: input forward wave A(0)=1, backward at output end B(L)=0
    % We will integrate from z=L down to 0 for B known at z=L, or from 0->L
    % A convenient approach: integrate coupled ODEs from z=L to 0 with boundary B(L)=0 and find A(0)
    % But simpler and common: integrate from z=0 to L with boundary A(0)=1 and assume B(L)=0
    % The second approach is approximate but acceptable for reflective grating strong enough.
    
    A = 1 + 0i;
    B = 0 + 0i;
    for iz = 1:Nz-1
        kz = kappa_z(iz);
        dz_local = dz;
        % RK4 step
        z_local = z(iz);
        % define derivative function
        f = @(Ain,Bin,zz,kappa,delta) deal( 1i*kappa*Bin.*exp(2i*delta*zz), ...
                                            -1i*kappa*Ain.*exp(-2i*delta*zz) );
        % k1
        [k1A,k1B] = f(A,B,z_local,kz,delta_z(iz));
        A2 = A + 0.5*dz_local*k1A;
        B2 = B + 0.5*dz_local*k1B;
        % k2 (use kappa at midpoint approx)
        kz2 = kappa_z(min(iz+1,Nz));
        [k2A,k2B] = f(A2,B2,z_local+0.5*dz_local,kz2,delta_z(min(iz+1,Nz)));
        A3 = A + 0.5*dz_local*k2A;
        B3 = B + 0.5*dz_local*k2B;
        % k3
        [k3A,k3B] = f(A3,B3,z_local+0.5*dz_local,kz2,delta_z(min(iz+1,Nz)));
        A4 = A + dz_local*k3A;
        B4 = B + dz_local*k3B;
        % k4
        kz4 = kappa_z(min(iz+2,Nz));
        [k4A,k4B] = f(A4,B4,z_local+dz_local,kz4,delta_z(min(iz+2,Nz)));
        % combine
        A = A + (dz_local/6)*(k1A + 2*k2A + 2*k3A + k4A);
        B = B + (dz_local/6)*(k1B + 2*k2B + 2*k3B + k4B);
    end
    % After integration, the backward wave at input (z=0) approximated by B
    R(iw) = abs(B)^2;  % reflected power (assuming input forward power 1)
end

% Plot results
figure('Name','CFBG Results','Position',[200 200 800 600]);
subplot(2,1,1);
plot((wvl_span*1e9), R, 'LineWidth',1.2);
xlabel('Wavelength (nm)'); ylabel('Reflectivity |R|^2');
title(['CFBG Reflectivity - Apodization: ' apod_type]);
grid on;

subplot(2,1,2);
plot(z*1e3, g, 'LineWidth',1.2);
xlabel('z (mm)'); ylabel('apodization g(z)');
title('Apodization profile along grating (normalized)');
grid on;

fprintf('Done. Peak reflectivity = %.3f at lambda = %.3f nm\n', max(R), wvl_span(R==max(R))*1e9);
