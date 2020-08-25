%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmitrii Briantcev, 2020
%
% Phase-screen generation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phz = Phz (params, time_step, C1, C2)

% Fried parameter for one step of propagation
r0_step = (0.423*params.k^2*params.Cn2*params.dz)^(-3/5); 

% Low frequency phase-screen generation
phz_l = phz_subharmonics(r0_step, params.N, params.delta, params.L0,...
   params.l0, C1, time_step, params.norm_velocity);

% High frequency phase-screen generation
phz_h = phz_high_freq(r0_step, params.N, params.delta, params.L0,...
   params.l0, C2, time_step, params.norm_velocity);

% Resulting 
phz = phz_l + phz_h;

end