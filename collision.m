%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collision.m: cillision step
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

A([8 9], (rho > rho_a))  = 1/tau_l;
A([8 9], (rho <= rho_a)) = 1/tau_v;

u_squ = ux.^2 + uy.^2;
me(1, :) = ones(1, lxy);
me(2, :) = -2 + 3*u_squ;
me(3, :) = 1 - 3*u_squ;
me(4, :) = ux;
me(5, :) = -ux;
me(6, :) = uy;
me(7, :) = -uy;
me(8, :) = ux.^2 - uy.^2;
me(9, :) = ux.*uy;
me       = ones(9, 1)*rho.*me;

m  = M*reshape(ff, 9, lxy);
ms = m - A.*(m - me) + (1 - A/2).*S;

fe = reshape(Minv*ms, 9, lx, ly);
% fe(:, obst) = ff(:, obst); % no collision on boundary nodes