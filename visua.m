%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visua.m: visualization
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

if (mod(time, Nwri) == 1)
%     save(['time', num2str(time)], 'fn');
    subplot(2, 2, 1);
    imagesc(flipud(reshape(rho, lx, ly)'));
    colorbar
    title('Density');
    axis equal off; drawnow
    subplot(2, 2, 2);
    quiver(x, y, reshape(ux, lx, ly), reshape(uy, lx, ly));
    colorbar
    title('Velocity');
    axis equal off; drawnow
    subplot(2, 2, 3);
    imagesc(flipud(reshape(p, lx, ly)'));
    colorbar
    title('Pressure');
    axis equal off; drawnow
    subplot(2, 2, 4);
    imagesc(flipud(reshape(T, lx, ly)'));
    colorbar
    title('Temperature');
    axis equal off; drawnow
end
