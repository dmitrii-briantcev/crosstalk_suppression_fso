%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Title: ft2
%    Author: Jason D. Schmidt.
%    Date: 2010
%    Availability: Numerical Simulation of Optical Wave Propagation
%                  With examples in MATLAB®
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = ft2(g, delta)

G = fftshift(fft2(fftshift(g)))*delta^2;

end