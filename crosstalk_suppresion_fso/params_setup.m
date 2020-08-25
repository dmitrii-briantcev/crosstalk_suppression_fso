%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmitrii Briantcev, 2020
%
% Define the parameters of the experiment, system, turbulence 
% and the simulator operations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = params_setup(r0, SNR_dB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External definitions of [SNR_dB, r0] pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.r0_vec =     [  0.03,   0.05,  0.07,  0.09,  0.11];
params.SNR_dB_vec = [   25,     25,    25,    25,    25 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.SNR_dB = SNR_dB; % Signal-to-noise ratio
params.SNR = 10.^(params.SNR_dB/10);
params.lambda = 1550e-9; % Wavelength
params.k = 2*pi/params.lambda; % Wavenumber
params.w0 = 0.025; % Beam waist in meters
params.z = 1000;% Propagation distance in meters 
params.dz = 10; % Prop. step in meters (must be a divisor of params.z)
params.z_div = params.z/params.dz; % Number of steps
params.delay = 1; % Delay in the number of time steps
params.time_lim = 250; % Number of time steps plus params.delay

% Turbulence parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.r0 = r0; % Fried parameter
params.L0 = 100; % Outer scale 
params.l0 = 0.01; % Inner scale
% Structure refractive index parameter 
params.Cn2 = (0.0598825*params.lambda^2)/(params.r0^(5/3)*params.z);
params.wind_velocity = 5; % Wind speed
params.refresh_rate = 500; % HZ - peak refresh rate of an SLM
% Wind velocity in [m/(SLM tick time)]
params.norm_velocity = params.wind_velocity/params.refresh_rate; 

% Reciver and grid parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.D = 0.5; % Length of one side of square phase screen [m]
% Number of grid points per side, should be a power of 2 for performance
params.N = 256; 
params.delta = params.D/params.N; % Distance between nodes
params.x = ((-params.N)/2:(params.N)/2-1)*params.delta; % x - grid
params.y = ((-params.N)/2:(params.N)/2-1)*params.delta; % y - grid
params.phase_error = 0.5; % Phase error of the reciever
params.responsivity = 0.8; % Electro-optical power conversion coefficient
 
% Simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.beamtype = 'LG'; % Can add more beamtypes (add as modes in OptMode.m)
params.mode = 1; % Decomposition mode. 1 - trapz, 2 - overlap sum.
params.iter = 1; % Number of independent runs

% Mode set definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  m and n should be of the same length
params.m_alph = [1,  0, 0, 0]; % Radial index
params.n_alph = [0, -1, 1, 0]; % Azimuthal index

end
