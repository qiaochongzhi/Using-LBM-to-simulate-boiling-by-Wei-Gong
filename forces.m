function [S, Fx, Fy] = forces( p, ff )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forces.m: calculate the interaction force between fluid nodes and the
%           interaction force between solid and fluid nodes
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

global lx ly lxy rho cc c_squ g tau_e tau_t ux uy ex ey sigm obst;

G = -1; 
w = single([4/3 1/3 1/3 1/3 1/3 1/12 1/12 1/12 1/12]);
S = single(zeros(9, lxy));
psx_re = zeros(9, lx, ly);

psx = sqrt(abs(2.*(p - rho.*c_squ)./G./cc./cc));
psx = reshape(psx, lx, ly);
for i = 1:9
       psx_re(i, :, :) = circshift(psx, [-ex(i), -ey(i)]);
end
psx = reshape(psx, 1, lxy);
psx_re = reshape(psx_re, 9, lxy);
Fmx = -G.*psx.*(w*(psx_re.*(ex'*ones(1, lxy))));
Fmy = -G.*psx.*(w*(psx_re.*(ey'*ones(1, lxy))));
Fmx(obst) = 0;
Fmy(obst) = 0;

Fbx = zeros(1, lxy);
Fby = (rho - mean(rho)).*(-g);

Fx = Fmx + Fbx;
Fy = Fmy + Fby;

ux     = (ex*reshape(ff, 9, lxy) + Fx/2)./rho;   
uy     = (ey*reshape(ff, 9, lxy) + Fy/2)./rho;
uF      = ux.*Fx + uy.*Fy;
Fm2     = Fmx.^2 + Fmy.^2; 
psx2    = psx.^2;
S(1, :) = 0;
S(2, :) = 6*uF + sigm*Fm2./psx2./(tau_e - 0.5);
S(3, :) = -6*uF - sigm*Fm2./psx2./(tau_t - 0.5);
S(4, :) = Fx;
S(5, :) = -Fx;
S(6, :) = Fy;
S(7, :) = -Fy;
S(8, :) = 2*(ux.*Fx - uy.*Fy);
S(9, :) = ux.*Fy + uy.*Fx;

end



