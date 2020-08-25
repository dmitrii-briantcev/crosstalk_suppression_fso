%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmitrii Briantcev, 2020
%
% Returns full decomposition vector of the disturbed signal
% E - disturbed field, CoeffArray - vector of decomposition coefficients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CoeffArray = Decompose (E, params)

m = params.m_alph;
n = params.n_alph;

% Array initialization
CoeffArray = zeros(length(m),1);

% Coefficient calculation
for i = 1:length(params.m_alph)
    E_decomp = clear_prop(m(i), n(i), params);
    CoeffArray(i) = Coeff(E,E_decomp,params);
end

end