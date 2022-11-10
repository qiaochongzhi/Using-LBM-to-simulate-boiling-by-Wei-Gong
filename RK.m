function [ T_new ] = RK( T_old, rho, ux, uy )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RK.m: calculate the temperature distribution using fourth-order
%       Runge-Kutta scheme
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

global lx ly lxy R b cv lambda obst; 

ux     = reshape(ux, lx, ly);
uy     = reshape(uy, lx, ly);
T_old  = reshape(T_old, lx, ly);
rho    = reshape(rho, lx, ly);

h1 = K(T_old);
h1(obst) = 0;  % boundary temperature keep constant
h2 = K(T_old + h1/2);
h2(obst) = 0;
h3 = K(T_old + h2/2);
h3(obst) = 0;
h4 = K(T_old + h3);
h4(obst) = 0;
T_new = reshape(T_old + (h1 + 2*h2 + 2*h3 + h4)/6, 1, lxy);

    function Kt = K(T)
        Kt = -(ux.*partx(T) + uy.*party(T)) + part2(T)./rho./cv - ...
            (partx(ux) + party(uy)).*T.*R./cv./(1 - b.*rho);
    end

    function px = partx(fx)
        % df/dx
        fxl = circshift(fx, [1, 0]);
        fxr = circshift(fx, [-1, 0]);
        px  = (fxr - fxl)/2;
    end

    function py = party(fy)
        % df/dy
        fyu = circshift(fy, [0, -1]);
        fyd = circshift(fy, [0, 1]);
        py  = (fyu - fyd)/2;
    end

    function p2 = part2(f)
        % d(lambda*df/dx)/dx + d(lambda*df/dy)/dy
        fxl = circshift(f, [1, 0]);
        fxr = circshift(f, [-1, 0]);
        fyu = circshift(f, [0, -1]);
        fyd = circshift(f, [0, 1]);
        lal = circshift(lambda, [1, 0]);
        lar = circshift(lambda, [-1, 0]);
        lau = circshift(lambda, [0, -1]);
        lad = circshift(lambda, [0, 1]);
        p2  = (lambda + lar)./2.*(fxr - f) - ...
            (lambda + lal)./2.*(f - fxl) + ...
            (lambda + lau)./2.*(fyu - f) - ...
            (lambda + lad)./2.*(f - fyd);
    end

end

