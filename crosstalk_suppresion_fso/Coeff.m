%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmitrii Briantcev, 2020
%
% Calculate generalized fourier coefficient (E, Ba, mode, params)
% E - field, Ba - basis function
% mode = 1 - trapezoidal intagration decompose
% mode = 2 - sum overlap decompose
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Coeff = Coeff(E, Ba, params)

if params.mode == 1  
    Coeff1 = trapz(params.y,trapz(params.x,E.*conj(Ba)), 2);
    Coeff2 = trapz(params.y,trapz(params.x,Ba.*conj(Ba)),2); 
elseif params.mode == 2  
    Coeff1 = sum(sum(E.*conj(Ba)), 2);
    Coeff2 = sum(sum(Ba.*conj(Ba)),2);  
end

Coeff = Coeff1/Coeff2;

end