clear variables;
% Numerical simulation as described in Mumford (2020)


% domain
Grid.x = 5;     Grid.dx = 0.1;      Grid.Nx = Grid.x/Grid.dx + 1;
Grid.z = 5;     Grid.dz = 0.05;     Grid.Nz = Grid.z/Grid.dz + 1;

x = linspace(0, Grid.x, Grid.Nx);
z = linspace(0, Grid.z, Grid.Nz);

t = 0;                  % start time
t_end = 10000;             % end time
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

Fluid.S_r = 0.13;       % residual wetting saturation
Fluid.S_gcr = 0.15;     % critical gas saturation

Fluid.sigma = 0.0623;   % interfacial tension of air/water

Fluid.P_cdim = 0.18557;                         % dimless cap pressure
Fluid.P_D = Fluid.P_cdim*Fluid.sigma*...
    (Fluid.por/Fluid.k)^(0.5);                  % displacement pressure

P_w = (Fluid.rho_w/1000)*9.8*Grid.z + 1.01*10^5;     % water pressure

% initial values
% initial temperature (Kelvin)
T = (10+273.15)*ones(Grid.Nz,Grid.Nx);     

% sink term that represents the heat comsumed by co-boiling
Q = zeros(Grid.Nz, Grid.Nx);      

% Initial saturations
S_g = zeros(Grid.Nz,Grid.Nx);         % initial gas saturation
S_n = 0.01*ones(Grid.Nz,Grid.Nx);     % initial water saturation
S_w = ones(Grid.Nz,Grid.Nx) - S_n;    % initial NAPL saturation
%%
% heat flux due to the heaters
f = zeros(Grid.Nz, Grid.Nx);
f_l = -(20/0.15);
f_r = (20/0.15);

% intital local gas pressure
P_g = ((S_w - Fluid.S_r*ones(size(S_w)))./((1-Fluid.S_r)*ones(size(S_w))...
    -S_w)).^(-1/Fluid.Lambda).*Fluid.P_D + P_w;

% initial K_e
K_e = (kappa*(ones(size(S_g)) - S_g))./(1 + (kappa - 1)*(1 - S_g));      

% initial thermal conductivity
lambda = K_e*(lambda_sat - lambda_dry) + lambda_dry;      

% initial heat capacity
heat_cap = S_w*Fluid.por*Fluid.rho_w*Fluid.C_pw + ...
    S_n*Fluid.por*Fluid.rho_n*Fluid.C_pn + ...
    (1-Fluid.por)*Fluid.rho_s*Fluid.C_ps; 

% volume of water and NAPL
V_w = S_w * Fluid.por * Grid.dx * Grid.dz;
V_n = S_n * Fluid.por * Grid.dx * Grid.dz;

% molar mass of water: 18.01528 g/mol
% molar mass of 1,1,1-Trichloroethane: 133.4 g/mol
% initial moles of water and NAPL
n_w = Fluid.rho_w*V_w / 18.01528;     % moles of water
n_n = Fluid.rho_n*V_n / 131.4;      % moles of NAPL

n_gn = zeros(Grid.Nz, Grid.Nx);
n_gw = zeros(Grid.Nz, Grid.Nx);

times = [t];
temps = [mean2(T - 273.15*ones(Grid.Nz,Grid.Nx))];

log_array = zeros(size(T));

%%
while t < 86400
    
    % compute temp
    T = temp_v3(Grid, T, Q,lambda, heat_cap,f_l, f_r);
    
    % compute vapor pressure using Antoine eqn
    P_nv = exp(19.796*ones(size(T)) - 2289.8*ones(size(T))...
        ./(T - 83.445*ones(size(T))));
    P_wv = exp(23.195*ones(size(T)) - 3814*ones(size(T))...
        ./(T - 46.29*ones(size(T))));
    
    % check if T reaches co-boiling temp
    if max(max(P_wv + P_nv)) >= (P_w + Fluid.P_D)
        
         % Compute the heat source/sink term
         % Issues with Q, so need to look over this again
         Q = gradient(lambda.* gradient(T));
        
         log_array = P_wv + P_nv >= (P_w + Fluid.P_D);
         
         Q = Q .*log_array;
         
         %P_wv = P_wv .* log_array; P_nv = P_nv .* log_array;
         
         n_gn = n_gn + Q./ (Fluid.L_w*(P_wv./P_nv) +...
             Fluid.L_n*ones(size(Q)));
%         n_gn = 0.0015*Grid.dt*n_gn;
         n_gw = n_gw + n_gn.*(P_wv./P_nv);
         
         n_gn = n_gn .* log_array; n_gw = n_gw .* log_array;
        
         
         % update moles and volume of NAPL and water (liquid phase)

         n_n = n_n - n_gn;
         n_n = n_n .* (n_n >= 0);
        
         n_w = n_w - n_gw;
         n_w = n_w .* (n_w >= 0);

         V_n = 131.4*n_n / Fluid.rho_n;
         V_w = 18.01528*n_w / Fluid.rho_w;
         
         % compute Q
%          Q = (Fluid.L_w * (n_gw+n_gn) + Fluid.L_n * (n_gw+n_gn)) /...
%              (Grid.dx * Grid.dz * Grid.dt);
         
         % compute volume of gas using the ideal gas law
         % NOTE: These are incorrect and need to be changed
         V_gn = 8.314462*(n_gn.*T) ./ P_nv;         % NAPL vapor
         V_gw = 8.314462*(n_gw.*T) ./ P_wv;         % water vapor
         
         % update saturations
         S_n = V_n / 0.0015;
         S_w = V_w / 0.0015;
         S_g = ones(size(S_g)) - (S_n + S_w);
%         S_g = (V_gn + V_gw)/0.0015;
         
         % macro-IP
         if max(max(S_g)) > Fluid.S_gcr
             %break
             clusters = propagate(S_g, Fluid.S_gcr);
         end
         
         % update heat capacity
         heat_cap = S_w*Fluid.por*Fluid.rho_w*Fluid.C_pw + ...
             S_n*Fluid.por*Fluid.rho_n*Fluid.C_pn + ...
             (1-Fluid.por)*Fluid.rho_s*Fluid.C_ps;

         % update thermal conductivity
         K_e = (kappa*(1 - S_g))./(1 + (kappa - 1)*(1 - S_g));    
         lambda = K_e*(lambda_sat - lambda_dry) + lambda_dry;

    end

    t = t + Grid.dt;
    
    times = [times; t];
    temps = [temps; mean2(T - 273.15*ones(Grid.Nz,Grid.Nx))];
end

T = T - 273.15*ones(Grid.Nz,Grid.Nx);

figure(1)
[xx,zz]=meshgrid(x,z);
[cs, hs] = contourf(xx,zz,flip(T,1));
set(hs,'EdgeColor','none')
colorbar
caxis([10 600])
colormap(jet)

figure(2)
plot(times, temps, 'Linewidth', 4)
xlabel('t (seconds)')
ylabel('Temperature (celsius)')
set(gca, 'Fontsize', 20)
title('Average temperature')

% figure(3)
% subplot(2,1,1);
% contourf(xx,zz,S_g,[.1 .2 .5]);
% 
% subplot(2,1,2);
% contourf(xx,zz,S_g,[.1 .2 .5]);
% 
% hold on;
% N = size(clusters);
% 
% for i = 1:N(2)
%     plot(xx(clusters(:,i)==1),zz(clusters(:,i)==1), '.', 'MarkerSize',20)
% end