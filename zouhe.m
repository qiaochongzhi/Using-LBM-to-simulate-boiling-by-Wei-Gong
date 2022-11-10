function [ ff ] = zouhe( ff )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zouhe.m: Zou-He boundary condition
%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shan Chen Lattice Boltzmann sample in Matlab
% Copyright Wei Gong
% Address: Nottingham NG7 2RD, UK
% E-mail: wei.gong@nottingham.ac.uk
% Reference: Li, Qing, et al. "Lattice Boltzmann modeling of boiling heat 
%            transfer: The boiling curve and the effects of wettability." 
%            International Journal of Heat and Mass Transfer 85 (2015): 
%            787-796.

global obst_b obst_u

ff(3, obst_b) = ff(5, obst_b);
ff(6, obst_b) = ff(8, obst_b) + (ff(4, obst_b) - ff(2, obst_b))/2;
ff(7, obst_b) = ff(9, obst_b) + (ff(2, obst_b) - ff(4, obst_b))/2;

ff(5, obst_u) = ff(3, obst_u);
ff(8, obst_u) = ff(6, obst_u) + (ff(2, obst_u) - ff(4, obst_u))/2;
ff(9, obst_u) = ff(7, obst_u) + (ff(4, obst_u) - ff(2, obst_u))/2;

end

