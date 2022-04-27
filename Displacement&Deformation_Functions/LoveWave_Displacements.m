function [R_out] = LoveWave_Displacements(r_in,S)
%
% Author: JessyJP (2022) % License: GPLv3 @ LICENCE.md
%

    % Unpack coordinates 
    x = r_in(:,1);
    y = r_in(:,2);
    z = -r_in(:,3);
    
    % Parameters INPUT
    c = S.c;
    c1 = S.c1;
    c2 = S.c2;
    A1 = S.A1;
    A2 = S.A2;
    h  = S.h;
    omega =S.omega;
    t = S.t;
    
    % Preallocate
    Uy = zeros(size(y));   
    % Displacements
    % Y translation
    if any(z<h)
    Uy(z<h)  = A1*cos( omega*z(z<h)*sqrt(1/c1^2 - 1/c^2)).*sin(omega*(t-x/c));
    end
    if any(h<=z)
    Uy(h<=z) = A2*exp(-omega*z(h<=z)*sqrt(1/c^2 - 1/c2^2)).*sin(omega*(t-x/c));
    end
        
        
    % X translation
    R_out(:,1) = x;
    % Y translation
    R_out(:,2) = y + (Uy);
    % Z translation
    R_out(:,3) = -z;
end