% extract the particle parameters
% Author: Tiberius Georgescu

clear all;
close all;
clc;

% Experimental constants
T0 = 273.15 + 20;   % lab air temperature, K, INSERT BY HAND
p_mbar = 6.7e0;     % pressure in mbar, INSERT BY HAND
p0 = p_mbar*100;    % pressure converted to Pa
%extracted params hpsd
file_name = 'C4 07_A1f0gn_out_ref_3.txt';   % name of file containing fit parameters,INSERT BY HAND
data = readmatrix(file_name);               % read the numbers from the .txt file
A1 = data(1);   % area under PSD
f0= data(2);    % central frequency
g = data(3);    % the damping rate
n = data(4);    % the noise floor
pkname = file_name(end-4);
omega_0 = 2*pi*f0;
Gamma_0 = 2*pi*g;
%% now extract the radius
rho_SiO2 = 1800;    % SiO2 particle, 1800 kg/m3, [1] 
M_m = 28.97e-3;     % molar mass of air
N_A = 6.022e23;     % Avogadro's number
k_B = 1.38e-23;     % Boltzmann's constant
m_gas = M_m/N_A;    % average mass of a molecule of gas
pre_factor = 0.619*9/...
            ( (2*pi)^0.5 * rho_SiO2 );
r_s =   pre_factor *...             % sphere radius
        ( M_m / (N_A*k_B*T0) )^0.5 * p0/Gamma_0;   
V_s = 4*pi*r_s.^3/3;    % sphere volume
m_s = rho_SiO2 * V_s;   % sphere mass
C1 = ( A1*m_s*omega_0^2/(k_B*T0) )^0.5; %conversion factor, V/m
%% write the extracted parameters
% T0 p0 rho(assumed) omega_0 Gamma_0 r_s m_s C1
extracted_params = [
    T0 % lab air temperature
    p0 % lab air pressure
    rho_SiO2 % assumed particle density
    A1 % area under the PSD 
    f0 % oscillation frequency of translational mode in Hz
    g  % damping factor in Hz
    n  % noise floor approximate level 
    omega_0     % frequency of translational mode in rad/s
    Gamma_0     % damping factor in rad/s
    r_s         % (spherical) particle radius
    m_s         % (spherical) particle mass
    C1          % conversion factor, V/m
    ];
dlmwrite([file_name(1:2),' ', file_name(21:22), '_TprhoA1f0gnO0G0rmC1_',pkname,'.txt'], extracted_params, 'precision', 15);
% the line above writes a .txt file containing the extracted parameters

% References:
% [1] The SiO2 particles I am using are from Corpuscular:
% https://microspheres-nanospheres.com/
% The density measured by the manufacturer is 1.8 g/cm3. 
% I obtained this value by requesting it from the manufacturer.
% One must be aware that Blakemore et al.[2] have measured the mass and
% density of three SiO2 particles of sizes: 2.348, 2.345 and 2.355 um
% and found their respective densities to be: 1.550, 1.554 and 1.562 g/cm3.
% Therefore, the 1.8 g/cm3 density assumption might not be accurate.
% [2] Blakemore, Charles P., et al. "Precision mass and density measurement of individual optically levitated microspheres." Physical Review Applied 12.2 (2019): 024037.
% [3] Hebestreit, Erik. Thermal properties of levitated nanoparticles. Diss. ETH Zurich, 2017.