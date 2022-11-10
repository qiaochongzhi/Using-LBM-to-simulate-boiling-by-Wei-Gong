%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constant.m: consant setting 
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

global lx ly lxy R b cv rho rho_a lambda cc c_squ tau_e tau_t ux uy ...
    ex ey sigm;

% D2Q9 lattice constants
lx     = single(200);        % x length of computational domain
ly     = single(150);        % y length of computational domain
lxy    = lx*ly;
cc     = single(1);          % lattice speeds
c_squ  = cc*cc/3;            % square of sound speed
ex     = single([0 1 0 -1 0 1 -1 -1 1]);                % velocity scheme
ey     = single([0 0 1 0 -1 1 1 -1 -1]);
M      = single([ 1  1  1  1  1  1  1  1  1;
                 -4 -1 -1 -1 -1  2  2  2  2;
                  4 -2 -2 -2 -2  1  1  1  1;
                  0  1  0 -1  0  1 -1 -1  1;
                  0 -2  0  2  0  1 -1 -1  1;
                  0  0  1  0 -1  1  1 -1 -1;
                  0  0 -2  0  2  1  1 -1 -1;
                  0  1 -1  1 -1  0  0  0  0;
                  0  0  0  0  0  1 -1  1 -1]);
Minv   = single([ 4 -4  4  0  0  0  0  0  0;   % inverse matrix
                  4 -1 -2  6 -6  0  0  9  0;
                  4 -1 -2  0  0  6 -6 -9  0;
                  4 -1 -2 -6  6  0  0  9  0;
                  4 -1 -2  0  0 -6  6 -9  0;
                  4  2  1  6  3  6  3  0  9;
                  4  2  1 -6 -3  6  3  0 -9;
                  4  2  1 -6 -3 -6 -3  0  9;
                  4  2  1  6  3 -6 -3  0 -9])/36;
sigm   = single(1.2);
       
% parameters in YUAN C-S EOS
a   = single(3/49);
b   = single(2/21);
R   = single(1);
ome = single(0.344);
Tc  = 0.0778/0.45724*a/(b*R);   % critical temperature
Ts  = 0.86*Tc;                  % initial temperature
dT  = single(0.0137);
Tb  = Ts + dT;

% General flow parameters
rho_l  = single(6.4989);        % density
rho_v  = single(0.3797);
rho_a  = (rho_l + rho_v)/2;     % average density
nu_l   = single(0.1);           % viscosity
nu_v   = single(0.5/3);
tau_l  = 3*nu_l + 0.5;          % relaxation time
tau_v  = 3*nu_v + 0.5;
tau_e  = single(1/1.1);
tau_t  = single(1/1.1);
A      = single([1 1/tau_e 1/tau_t 1 1.1 1 1.1 1 1]'*ones(1, lxy));   % orthogonal transformation matrix
cv     = single(6);                     % heat capacity
rho    = single(zeros(1, lxy));
p      = single(zeros(1, lxy));         % pressure
phi    = single(zeros(1, lxy));
ux     = single(zeros(1, lxy));         % velocity
uy     = single(zeros(1, lxy));
T      = single(zeros(1, lxy));         % temperature 
Fx     = single(zeros(1, lxy));
Fy     = single(zeros(1, lxy));
lambda = single(zeros(lx, ly));         % heat conductivity
m      = single(zeros(9, lxy)); 
me     = single(zeros(9, lxy));
ff     = single(zeros(9, lx, ly));      % distribution function
fe     = single(zeros(9, lx, ly));      % distribution function







