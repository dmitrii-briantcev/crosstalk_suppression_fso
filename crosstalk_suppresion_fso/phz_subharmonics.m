%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmitrii Briantcev, 2020
%
% Generation of the low frequency correction for Fourier-generated
% phase-screens with time dependence
%
% Based on:
%    Title: ft_sh_phase_screen
%    Author: Jason D. Schmidt.
%    Date: 2010
%    Availability: Numerical Simulation of Optical Wave Propagation
%                  With examples in MATLAB®
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phz_low = phz_subharmonics(r0, N, delta, L0, l0,...
    C1, time_step, velocity)

% Shifted spatial grid
% sqrt(2) for diagonal movement over the screen, provides longer
% propagation without phasescreen repetition
[x, y] = meshgrid((-N/2 : N/2-1) * delta + velocity*time_step/sqrt(2));

% Low-frequency screen compensation initialization
phz_low = zeros(N);

for p = 1:3
    
    df = 1/(3^p*N*delta); % Frequency grid spacing 
    [fx, fy] = meshgrid((-1:1)*df); % Frequency grid 
    [~, f] = cart2pol(fx, fy); % Polar grid
    
    % Modified von Karman atmospheric phase PSD
    VK = 0.023*r0^(-5/3)*exp(-(f/(5.92/l0/(2*pi))).^2)./(f.^2+(1/L0)^2).^(11/6);
    VK(2,2) = 0; % Remove constant term
    
    % Random draws of Fourier coefficients
    RD = C1{p}.*sqrt(VK)*df;
    
    % SH initialization
    SH = zeros(N);
    
    % Accumulate SH
    for j = 1:9
        SH = SH + RD(j)*exp(1i*2*pi*(fx(j)*x+fy(j)*y));
    end
    
    % Accumulate phase-screen
    phz_low = phz_low + SH;
    
end

% Center the lo-freq shift
phz_low = real(phz_low) - mean(real(phz_low(:)));

end