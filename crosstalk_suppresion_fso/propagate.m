%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmitrii Briantcev, 2020
%
% Propagation of the arbitrary optical field to a distance dz
%
% Based on:
%    Title: ang_spec_prop
%    Author: Jason D. Schmidt.
%    Date: 2010
%    Availability: Numerical Simulation of Optical Wave Propagation
%                  With examples in MATLAB®
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E_prop = propagate(E,params,dz)

% Generate frequency grid 
N = params.N;
df = 1/(params.N*params.delta);
[fx, fy] = meshgrid((-N/2:N/2-1)*df);
fr2 = fx.^2 + fy.^2;

% Quadratic phase factor
Q = exp(-1i*pi^2*2*dz/params.k*fr2 + 1i*params.k*dz);

% Compute the propagated field
E_prop = ift2(Q.*ft2(E,params.delta),df);

end






