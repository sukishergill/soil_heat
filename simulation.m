clear variables;
% Numerical simulation as described in Mumford (2020)

% domain
Grid.x = 5;     Grid.dx = 0.1;       Grid.Nx = Grid.x/Grid.dx + 1;
Grid.z = 5;     Grid.dz = 0.05;      Grid.Nz = Grid.z/Grid.dz + 1;

x = linspace(0, Grid.x, Grid.Nx);
z = linspace(0, Grid.z, Grid.Nz);

[xx,zz]=meshgrid(x,z);

t = 0;                     % start time
t_end = 10000;             % end time
Grid.dt = 720;             % time step (seconds)

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

% Fluid.S_r = 0.13;       % residual wetting saturation
% Fluid.S_gcr = 0.15;     % critical gas saturation

Fluid.sigma = 0.0623;   % interfacial tension of air/water

P_cdim = 0.18557;                                   % dimless cap pressure
Fluid.P_D = P_cdim*Fluid.sigma*...
    (Fluid.por/Fluid.k)^(0.5);                      % displacement pressure

P_w = 9.8*(Fluid.rho_w/1000)*(ones(Grid.Nz,Grid.Nx)...
    .*z') + 1.01*10^5*ones(Grid.Nz,Grid.Nx);          % water pressure

% initial values
% initial temperature (Kelvin)
T = (10+273.15)*ones(Grid.Nz,Grid.Nx);     

% sink term that represents the heat comsumed by co-boiling
Q = zeros(Grid.Nz, Grid.Nx);    

% Initial saturations
S_g = zeros(Grid.Nz,Grid.Nx);           % initial gas saturation
S_n = 0.01*ones(Grid.Nz,Grid.Nx);       % initial water saturation
S_w = ones(Grid.Nz,Grid.Nx) - S_n;      % initial NAPL saturation
Fluid.S_r = 0.13;                             % residual wetting saturation
Fluid.S_gcr = 0.15;                           % critical gas saturation

% heat flux due to the heaters
f_l = -10/(0.15);
f_r = 0;

% initial K_e
K_e = (kappa*(ones(size(S_g)) - S_g))./(1 + (kappa - 1)*(1 - S_g));   

%ep = rand(size(T));
ep = 0 * ones(size(T));

% initial thermal conductivity
lambda = K_e*(lambda_sat - lambda_dry) + lambda_dry + ep;   

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

times = [t];
temps = [mean2(T - 273.15*ones(Grid.Nz,Grid.Nx))];

%%
while t < 500000
    
    % compute temp
    T = temp_v3(Grid, T, Q,lambda, heat_cap,f_l, f_r);
    
    % compute vapor pressure using Antoine eqn
    P_nv = exp(19.796*ones(size(T)) - 2289.8*ones(size(T))...
        ./(T - 83.445*ones(size(T))));
    P_wv = exp(23.195*ones(size(T)) - 3814*ones(size(T))...
        ./(T - 46.29*ones(size(T))));
    
    co_boil = (P_wv + P_nv) >= ((P_w + Fluid.P_D*ones(size(T))) + ...
        10*cos(t/100)*ones(size(T)));
    
    % check if T reaches co-boiling temp
    if any(co_boil, 'all') == 1
        
         % Compute the heat source/sink term
         Q = gradient(lambda.* gradient(T)) .* co_boil;

         % compute the moles of vapor produced using the energy
         % balance equation and Dalton's law
         n_gn = (Q./ (Fluid.L_w*(P_wv./P_nv) +...
             Fluid.L_n*ones(size(Q))));                 % energy balance
         n_gw = n_gn.*(P_wv./P_nv);                     % Dalton's law
         
         
         % Compute capillary and gas pressures
         % capillary pressure
         P_c = Fluid.P_D * ((max(S_w,Fluid.S_r) - Fluid.S_r*...
             ones(size(S_w)))./((1-Fluid.S_r)...
             *ones(size(S_w))-S_n)).^(-1/Fluid.Lambda);

         % local gas pressure   
         P_g =  P_c + P_w;

         % compute volume of gas using the ideal gas law
         V_gn = (8.314462*(n_gn.*T) ./ P_g);       % NAPL vapor
         
         if min(min(V_n - V_gn)) < 0
            V_gn =  V_gn .* (V_n > V_gn) + V_n .* (V_n < V_gn);
         end
         
         V_gw = (8.314462*(n_gw.*T) ./ P_g);       % water vapor
         
         if min(min(V_w - V_gw)) < 0
            V_gw =  V_gw .* (V_w > V_gw) + V_w .* (V_w < V_gw);
         end
         
         V_n = (V_n - V_gn) .* (V_n >= V_gn);
         V_w = (V_w - V_gw) .* (V_w >= V_gw);
      
         % update saturations
         S_n = V_n / 0.0015;
         S_w = V_w / 0.0015;
         S_g = S_g + (V_gw + V_gn)/0.0015;
              
         % macro-IP
%          if max(max(S_g)) > Fluid.S_gcr
%              clusters = propagate(S_g, Fluid.S_gcr);
%              break
%          end
         
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