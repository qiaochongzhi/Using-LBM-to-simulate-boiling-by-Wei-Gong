%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boiling.m: the main program for simulating droplet motion on surfaces
%            with Shan Chen lattice Boltzmann model
%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shan Chen Lattice Boltzmann sample in Matlab
% Reference: Li, Qing, et al. "Lattice Boltzmann modeling of boiling heat 
%            transfer: The boiling curve and the effects of wettability." 
%            International Journal of Heat and Mass Transfer 85 (2015): 
%            787-796.

clear, clc

global g;

format long;

t_max = single(1000000);    % maximum iteration
Nwri  = single(100);        % output data frequency

constant;         % constant setting
initialization;   % initial condition setting

% load time1;

% iteration (first 100 iters, without gravity)
g = single(0);      
for time = single(1:100000)    
    disp(time);
    % streaming
    for i = 1:9
        ff(i, :, :) = ...  
            circshift(fe(i, :, :), [0, ex(i), ey(i)]);
    end
    % Zou-He boundary conditions
    ff = zouhe(ff); 
    % macro parameters
    macrop;
    % interaction forces
    [S, Fx, Fy] = forces( p, ff ); 
    % collision
    collision; 
    % visualization
    visua;
end

% iteration (after 100 iters, with gravity)













