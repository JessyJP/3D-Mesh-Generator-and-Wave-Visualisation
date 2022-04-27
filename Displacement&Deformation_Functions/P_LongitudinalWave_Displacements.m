function [R_out] = P_LongitudinalWave_Displacements(r_in,S)
%
% Author: JessyJP (2022) % License: GPLv3 @ LICENCE.md
%

    % Unpack coordinates 
    x = r_in(:,1);
    y = r_in(:,2);
    z = r_in(:,3);
    
    % Parameters INPUT
    t = S.t;
    A = S.A;    
    k = S.k;
    omega =S.omega;
 
    % Displacements
    % X translation
    Ux = A*exp(1i*(omega*t-k*x));        
        
    % X translation
    R_out(:,1) = x + real(Ux);
    % Y translation
    R_out(:,2) = y;
    % Z translation
    R_out(:,3) = z;
end