function [R_out] = RayleighWave_Displacements(r_in,S)
%
% Author: JessyJP (2022) % License: GPLv3 @ LICENCE.md
%

    % Unpack coordinates 
    x = r_in(:,1);
    y = r_in(:,2);
    z = r_in(:,3);
    
    % Parameters INPUT
    c = S.c;
    cl = S.cl;
    ct = S.ct;
    A = S.A1;
    omega =S.omega;
    t = S.t;
    
    % Parameter Simplification
    s1 = -sqrt(c^2 - cl^2);
    s2 = -sqrt(c^2 - ct^2);   
    % Rayleigh wave Displacements
    % X translation
    Ux= A*c*(exp(-s1*z)-(2*s2*s1/(ct^2))*exp(-s2.*z)).*exp(1i*(c.*x-omega.*t-pi/2));
    % Y translation
    Uz= A*s1*(exp(-s1*z)-((2*c^2)/((2*c^2)-(ct^2)))*exp(-s2.*z)).*exp(1i.*(c.*x-omega.*t));
    
     
    
    % X translation
    R_out(:,1) = real(Ux)+x;
    % Y translation
    R_out(:,2) = y;
    % Z translation
    R_out(:,3) = real(Uz)+z;
end