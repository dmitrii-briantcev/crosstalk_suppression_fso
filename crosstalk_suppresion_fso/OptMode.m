%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmitrii Briantcev, 2020
%
% Choice of optical mode
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = OptMode (m, n, params)

% Here code can be modified to propagate any SL field,
% By stating an appropriate beam type at params.beamtype in params_setup,
% and adding the generating function (like lg() in this case) to the code folder

if params.beamtype == 'LG'
    E = lg(m, n, params);
else
    return 
end
