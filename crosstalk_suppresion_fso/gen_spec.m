%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmitrii Briantcev, 2020
%
% Generate spectrum random draws for turbulence over a period of time defined in
% params.time_lim
% Can be modified for time-varying spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C1, C2] = gen_spec(params,seed)

rng(seed) % Initialize rnd
screen_num_lim = params.z_div;
C1 = cell(screen_num_lim);
C2 = cell(screen_num_lim);

for j = 1:screen_num_lim
    
	c1_init{1} = (randn(3) + 1i*randn(3));
    c1_init{2} = (randn(3) + 1i*randn(3));
    c1_init{3} = (randn(3) + 1i*randn(3));
    c2_init = (randn(params.N) + 1i*randn(params.N));
    
    % Here spectrum can be modified ("Bubbling effect")
    C1{j} = c1_init;
    C2{j} = c2_init;
 
end

end