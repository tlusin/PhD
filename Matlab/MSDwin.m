function [win_MSD, MSD_coeffs] = MSDwin(N,L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametri funkcije
%--------------------------------------------------------------------------
% N.............število vzorènih toèk
% L.............red kosinusne okenske funkcije
%               L=0 - pravokotno okno
%               L=1 - Hannovo okno

% win_MSD.......toèke kosinusne okenske funckije
% MSD_coeffs....koeficienti kosinusne okenske funkcije
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    L=1;
end

coeff = zeros(1,L);

switch (L)
    case 0
        MSD_coeffs = 1;
    
    otherwise
        C = @(m, p) factorial(m)/(factorial(m-p)*factorial(p));
        a0 = C(2*L,L)/2^(2*L);
        ah = @(h) C(2*L,L-h)/2^(2*L-1);

        for h = 1:L
            coeff(h) = (-1)^h*ah(h);
        end
        
        MSD_coeffs = [a0 coeff];
end

x = (0:N-1)'*2.0*pi/N;
win_MSD = cos(x*(0:L))*MSD_coeffs';

end % function
