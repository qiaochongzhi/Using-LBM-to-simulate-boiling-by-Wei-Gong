%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% macrop.m: calculate all the macro parameters
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

T      = RK(T, rho, ux, uy);
rho    = reshape(sum(ff), 1, lxy);   
phi    = (1 + (0.37464 + 1.54226*ome - 0.26992*ome^2)*(1 - sqrt(T/Tc))).^2;
p      = rho.*R.*T./(1 - b.*rho) - ...
       a.*phi.*rho.^2./(1 + 2.*b.*rho - b.^2.*rho.^2);
lambda = 0.028*reshape(rho, lx, ly);
