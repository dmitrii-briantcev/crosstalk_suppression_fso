%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmitrii Briantcev, 2020
%
% Generation of the Fourier-generated phase-screens with time dependence
%
% Based on:
%    Title: ft_phase_screen
%    Author: Jason D. Schmidt.
%    Date: 2010
%    Availability: Numerical Simulation of Optical Wave Propagation
%                  With examples in MATLAB®
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phz_high = phz_high_freq(r0, N, delta, L0, l0,...
    C2, time_step, velocity)

% Generate frequency grid
[fx, fy] = meshgrid((-N/2:N/2-1)*(1/(N*delta)));
[~, f] = cart2pol(fx, fy); 

% Modified von Karman atmospheric phase PSD
VK = 0.023*r0^(-5/3)*exp(-(f/(5.92/l0/(2*pi))).^2)./(f.^2+(1/L0)^2).^(11/6);

% Remove constant term
VK(N/2+1,N/2+1) = 0;

% Shifted hi-freq spctrum
% sqrt(2) for diagonal movement over the screen, provides longer
% propagation without phasescreen repetition
cn = C2.*sqrt(VK)*(1/(N*delta)).*exp(1i*2*pi*velocity*time_step/sqrt(2).*(fx + fy));

% Final Phase-screen generation
phz_high = real(ift2(cn,1));

end