% simulation for cylindrical domain
clear variables;

% domain
Grid.r = 5;     Grid.dr = 0.1;      Grid.Nr = Grid.r/Grid.dr;

r = linspace(Grid.dr, Grid.r, Grid.Nr);
% center shell is not included (to avoid singularity at r=0)
%r = r(2:end);

t = 0;                  % start time
t_end = 10000;          % end time
Grid.dt = 720;          % time step (seconds)

% parameters
Fluid.por = 0.3;              % porosity
Fluid.k = 1.03151e-12;        % permeability

Fluid.C_pw = 4.184/1000;           % heat capacity of water
Fluid.C_pn = 0.958/1000;           % heat capacity of TCE
Fluid.C_ps = 0.8/1000;             % heat capacity of soil

Fluid.rho_w = 1000000*1;              % density of water
Fluid.rho_n = 1000000*1.46;           % density of TCE
Fluid.rho_s = 1000000*2.7;            % grain density

Fluid.L_w = 41.47;            % latent heat of vaporization of water
Fluid.L_n = 31.24;            % latent heat of vaporization of TCE

lambda_sat = 2.75/1000;      % thermal conductivity of saturated soil
lambda_dry = 0.15/1000;      % thermal conductivity of dry soil

kappa = 1.9;            % soil texture dependent parameter

Fluid.Lambda = 2.5;           % por size distribution index

Fluid.sigma = 0.0623;   % interfacial tension of air/water

P_cdim = 0.18557;                                   % dimless cap pressure
Fluid.P_D = P_cdim*Fluid.sigma*...
    (Fluid.por/Fluid.k)^(0.5);                      % displacement pressure

P_w = (9.8*(Fluid.rho_w/1000)+ 1.01*10^5)*...
    (ones(Grid.Nr,1));                               % water pressure

% initial values
% initial temperature (Kelvin)
T = (10+273.15)*ones(Grid.Nr,1);

% sink term that represents the heat comsumed by co-boiling
Q = zeros(Grid.Nr,1);

% Initial saturations
S_g = zeros(Grid.Nr,1);                     % gas saturation
S_n = 0.01 * ones(Grid.Nr,1);               % NAPL saturation
S_w = ones(Grid.Nr,1) - S_n;                % water saturation

K_e = (kappa*(ones(size(S_g)) - S_g))./ ...
    (1 + (kappa - 1)*(1 - S_g));

% thermal conductivity
lambda = K_e*(lambda_sat - lambda_dry) +...
    lambda_dry;                              

% heat capacity
heat_cap = S_w*Fluid.por*Fluid.rho_w*Fluid.C_pw + ...
    S_n*Fluid.por*Fluid.rho_n*Fluid.C_pn + ...
    (1-Fluid.por)*Fluid.rho_s*Fluid.C_ps;           

% heat flux due to heaters
f_R = -40/0.15;
f_inf = 0;

%%
while t < 500000
   
   % compute temp
   T = temp_rad_1D(Grid, T, Q, lambda, heat_cap, f_R, f_inf);
   
   t = t + Grid.dt;
   
end

T = T - 273.15;

th = linspace(0, pi/2, 100);
[rr, tt] = meshgrid(r, th);
T = repmat(T, 100, 1);
T = reshape(T, 50, 100);
T = T';

[x, y, v] = pol2cart(tt, rr, T);

figure(1)
h = surf(x, y, v)
set(h, 'edgecolor', 'none')
colorbar;
axis square

% h = polar(xx, yy);
% hold on
% contourf(xx, yy, T')