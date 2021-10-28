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

Fluid.P_cdim = 0.18557;                         % dimless cap pressure
Fluid.P_D = 5.31*10^2;                          % displacement pressure
P_w = Fluid.rho_w*9.8*Grid.z + 1.01*10^5;       % water pressure

% domain
Grid.x = 20;    Grid.dx = 0.1;      Grid.Nx = Grid.x/Grid.dx;
Grid.z = 5;     Grid.dz = 0.05;     Grid.Nz = Grid.z/Grid.dz;

t = 0;                  % start time
t_end = 1296000;        % end time
Grid.dt = 0.1;               % time step

% initial values
% initial temperature
T = 10*ones(Grid.Nz,Grid.Nx);     

% sink term that represents the heat comsumed by co-boiling
Q = zeros(Grid.Nz, Grid.Nx);      

% Initial saturations
S_g = zeros(Grid.Nz,Grid.Nx);         % initial gas saturation
S_w = 0.5*ones(Grid.Nz,Grid.Nx);      % initial water saturation
S_n = 0.5*ones(Grid.Nz,Grid.Nx);      % initial NAPL saturation

% initial moles of water and NAPL
n_w = zeros(Grid.Nz,Grid.Nx);     % moles of water
n_n = zeros(Grid.Nz,Grid.Nx);     % moles of NAPL

% heat flux due to the heaters
f = zeros(Grid.Nz, Grid.Nx);

f(:,1) = -20*ones(1,Grid.Nz);
f(:,end) = -1*f(:,1);
f(:,50) = 40*ones(1,Grid.Nz);
f(:,100) = f(:,50);
f(:,150) = f(:,50);


% initial K_e
K_e = (kappa*(1 - S_g))./(1 + (kappa - 1)*(1 - S_g));      

% initial thermal conductivity
lambda = K_e*(lambda_sat - lambda_dry) + lambda_dry;      

% initial heat capacity
heat_cap = S_w*Fluid.por*Fluid.rho_w*Fluid.C_pw + ...
    S_n*Fluid.por*Fluid.rho_n*Fluid.C_pn + ...
    (1-Fluid.por)*Fluid.rho_s*Fluid.C_ps;    


%%
while t <= t_end
    
    % compute temp
    T = temp(Grid, T, Q, lambda, heat_cap, f);
    
    % compute vapor pressure using Antoine eqn
    P_wv = 10.^(8.07131*ones(size(T)) - 1730.63*ones(size(T))...
        ./(233.426*ones(size(T)) + T));
    P_nv = 10.^(7.11*ones(size(T)) - 1367.05*ones(size(T))...
        ./(235.809*ones(size(T))+T));
    
    % check if T reaches co-boiling temp
    if max(P_wv + P_nv) >= (P_w + Fluid.P_D)
         % compute mole fractions
         y_n = P_nv ./ (P_nv + P_wv);       % NAPL mole fraction
         y_w = P_nw ./ (P_nv + P_wv);       % water mole fraction
         
         % compute moles of water vapor and NAPL steam produced
         n_gn = y_n .* n_n;     % NAPL vapor
         n_gw = y_w .* n_w;     % water vapor
         
         % update moles of NAPL and water (liquid phase)
         n_n = n_n - n_gn;
         n_w = n_w - n_gw;
         
         % compute Q
         Q = (Fluid.L_w * n_gw + Fluid.L_n * n_gn) /...
             (Grid.dx * Grid.dz * Grid.dt);
         
         % compute volume of gas using the ideal gas law
         V_n = (n_gn.*T) ./ P_nv;         % NAPL vapor
         V_w = (n_gw.*T) ./ P_wv;         % water vapor
         
         
         % !!!!!!!!!!!!!!!!!!!
         % INCOMPLETE
         % !!!!!!!!!!!!!!!!!!!
         % Need to express volume of gas in terms of saturation, which will
         % update S_g. Also need to update S_n and S_w using the moles of
         % gas produced and/or the equilibium equations S_n + S_w + S_g = 1        
         
         % update heat capacity
         heat_cap = S_w*Fluid.por*Fluid.rho_w*Fluid.C_pw + ...
             S_n*Fluid.por*Fluid.rho_n*Fluid.C_pn + ...
             (1-Fluid.por)*Fluid.rho_s*Fluid.C_ps;

         % update thermal conductivity
         K_e = (kappa*(1 - S_g))/(1 + (kappa - 1)*(1 - S_g));    
         lambda = K_e*(lambda_sat - lambda_dry) + lambda_dry;
         
    end

    t = t + Grid.dt;
end
