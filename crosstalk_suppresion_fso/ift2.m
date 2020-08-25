%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Title: ift2
%    Author: Jason D. Schmidt.
%    Date: 2010
%    Availability: Numerical Simulation of Optical Wave Propagation
%                  With examples in MATLAB®
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = ift2(G, delta_f)

N = size(G, 1);
g = ifftshift(ifft2(ifftshift(G)))*(N*delta_f)^2;

end