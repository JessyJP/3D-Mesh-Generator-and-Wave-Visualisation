function [R_out] = S_TransverseWave_Displacements(r_in,S)
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
    omega =S.omega;
    k = S.k;
 
    % Displacements
    % Y translation
    Uy = A*cos(omega*t-k*z); 
    % Z translation
    Uz = A*sin(omega*t-k*x); 
        
    % X translation
    R_out(:,1) = x;
    % Y translation
    R_out(:,2) = y + Uy*0;
    % Z translation
    R_out(:,3) = z + Uz*1; 
end