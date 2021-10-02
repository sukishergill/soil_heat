clear all; clc;
% Numerical simulation as described in Mumford (2020)

% parameters
Fluid.por = 0.3;                    % porosity

Fluid.C_pw = 4.184;           % heat capacity of water
Fluid.C_pn = 0.958;           % heat capacity of TCE
Fluid.C_ps = 0.8;             % heat capacity of soil

Fluid.rho_w = 1;              % density of water
Fluid.rho_n = 1.46;           % density of TCE
Fluid.rho_s = 2.7;            % grain density

Fluid.L_w = 41.47;            % latent heat of vaporization of water
Fluid.L_n = 31.24;            % latent heat of vaporization of TCE

lambda_sat = 2.75;      % thermal conductivity of saturated soil
lambda_dry = 0.15;      % thermal conductivity of dry soil

kappa = 1.9;            % soil texture dependent parameter

Fluid.Lambda = 2.5;           % por size distribution index

Fluid.S_r = 0.13;       % residual wetting saturation
Fluid.S_gcr = 0.15;     % critical gas saturation

Fluid.sigma = 0.0623;   % interfacial tension of air/water

Fluid.P_cdim = 0.18557;       % dimensionless capillary pressure
Fluid.P_D = 5.31*10^2;        % displacement pressure

% domain
Grid.x = 20;    Grid.dx = 0.1;      Grid.Nx = Grid.x/Grid.dx;
Grid.z = 5;     Grid.dz = 0.05;     Grid.Nz = Grid.z/Grid.dz;

Grid.t = 0;                  % start time
Grid.t_end = 1296000;        % end time
Grid.dt = 0.1;               % time step

% initial values
% initial temperature
T = 10*ones(Grid.Nz,Grid.Nx);     

% sink term that represents the heat comsumed by co-boiling
Q = zeros(Grid.Nz, Grid.Nx);      

% !!!!!!!!!!!!!!!!!!!
% INCOMPLETE
% !!!!!!!!!!!!!!!!!!!
% Initial water and NAPL saturations are not given
S_g = zeros(Grid.Nz,Grid.Nx);     % initial gas saturation
S_w = 0.5*ones(Grid.Nz,Grid.Nx);      % initial water saturation
S_n = 0.5*ones(Grid.Nz,Grid.Nx);      % initial NAPL saturation

% !!!!!!!!!!!!!!!!!!!
% INCOMPLETE
% !!!!!!!!!!!!!!!!!!!
% heat flux due to the heaters
f_ext = 0;        % exterior heaters
f_int = 0;        % interior heaters

% initial K_e
K_e = (kappa*(1 - S_g))./(1 + (kappa - 1)*(1 - S_g));      

% initial thermal conductivity
lambda = K_e*(lambda_sat - lambda_dry) + lambda_dry;      

% initial heat capacity
heat_cap = S_w*Fluid.por*Fluid.rho_w*Fluid.C_pw + ...
    S_n*Fluid.por*Fluid.rho_n*Fluid.C_pn + ...
    (1-Fluid.por)*Fluid.rho_s*Fluid.C_ps;                 

%% 
T = temp(Grid, T, Q, lambda, heat_cap, f_ext, f_int);

%%
[n_gn, n_gw] = NAPL_mole(Fluid, T);

%%
Sg = macroIP(S_w, S_g, S_n, Fluid, Grid, T);

%%
while t <= t_end
    
    % compute temp
    T = temp(Grid, T, Q, lambda, heat_cap, f_ext, f_int);
    
    % compute moles of NAPL vapor
    %[n_gn, n_gw] = NAPL_mole(Fluid, T);
    
    S_g = macroIP(S_w, S_g, S_n, Fluid, Grid, T);

    % !!!!!!!!!!!!!!!!!!!
    % INCOMPLETE
    % !!!!!!!!!!!!!!!!!!!
    % NEED TO UPDATE HEAT CAPACITY, S_n and S_w
    % use equilibrium equations S_n + S_w + S_g = 1
    
    % update thermal conductivity
    K_e = (kappa*(1 - S_g))/(1 + (kappa - 1)*(1 - S_g));    
    lambda = K_e*(lambda_sat - lambda_dry) + lambda_dry;

    t = t + dt;
end
