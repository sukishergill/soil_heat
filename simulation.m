clear variables;
% Numerical simulation as described in Mumford (2020)

% domain
Grid.z = 5;     Grid.dz = 0.05;     Grid.Nz = Grid.z/Grid.dz + 1;
Grid.x = 5;     Grid.dx = 0.1;       Grid.Nx = Grid.x/Grid.dx + 1;

x = linspace(0, Grid.x, Grid.Nx);
z = linspace(0, Grid.z, Grid.Nz);

[xx,zz]=meshgrid(x,z);

t = 0;                     % start time
t_end = 10000;             % end time
Grid.dt = 720;             % time step (seconds)

% parameters
Fluid.por = 0.3;              % porosity
Fluid.k = 1.03151e-12;        % permeability
% ln(k) has mean -27.6 and variance 2 (moderate heterogeneity)
% Fluid.k = exp(-27.6 + sqrt(2)*randn(Grid.Nz, Grid.Nx));

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

Fluid.Lambda = 2.5;           % pore size distribution index

Fluid.sigma = 0.0623;   % interfacial tension of air/water

P_cdim = 0.18557;                                   % dimless cap pressure
Fluid.P_D = P_cdim*Fluid.sigma*...
    (Fluid.por./Fluid.k).^(0.5);                    % displacement pressure

P_w = 9.8*(Fluid.rho_w/1000)*(ones(Grid.Nz,Grid.Nx)...
    .*z') + 1.01*10^5*ones(Grid.Nz,Grid.Nx);          % water pressure

% P_w = 9.8*(25000/1000)*(ones(Grid.Nz,Grid.Nx)...
%     .*z') + 1.01*10^5*ones(Grid.Nz,Grid.Nx); 

% P_w = 9.8*(25000/1000)*(ones(Grid.Nz,Grid.Nx)...
%      .*z' - 2.5) + 1.255*10^5*ones(Grid.Nz,Grid.Nx);

% uniform water pressure
% P_w = 9.8*(Fluid.rho_w/1000)*(ones(Grid.Nz,Grid.Nx))...
%     + 1.01*10^5*ones(Grid.Nz,Grid.Nx);                  % water pressure

V_cell = Fluid.por * Grid.dx * Grid.dz;

% initial values
% initial temperature (Kelvin)
T = (10+273.15)*ones(Grid.Nz,Grid.Nx);     

% sink term that represents the heat comsumed by co-boiling
Q = zeros(Grid.Nz, Grid.Nx);    
Q_old = Q;

% Initial saturations
S_g = zeros(Grid.Nz,Grid.Nx);           % initial gas saturation
S_n = 0.1*ones(Grid.Nz,Grid.Nx);       % initial water saturation
S_w = ones(Grid.Nz,Grid.Nx) - S_n;      % initial NAPL saturation
Fluid.S_r = 0.13;                             % residual wetting saturation
Fluid.S_gcr = 0.15;                           % critical gas saturation

% heat flux due to the heaters
% f_l = -80/2.75;
f_l = -10 / 0.15;
% f_r = 10/0.15;
% f_l = 0; 
f_r = 0;

% initial K_e
K_e = (kappa*(ones(size(S_g)) - S_g))./(1 + (kappa - 1)*(1 - S_g));   

% initial thermal conductivity
lambda = K_e*(lambda_sat - lambda_dry) + lambda_dry;   

% initial heat capacity
heat_cap = S_w*Fluid.por*Fluid.rho_w*Fluid.C_pw + ...
    S_n*Fluid.por*Fluid.rho_n*Fluid.C_pn + ...
    (1-Fluid.por)*Fluid.rho_s*Fluid.C_ps; 

% volume of water and NAPL
V_w = (S_w - Fluid.S_r) * V_cell;
V_n = S_n * V_cell;

% molar mass of water: 18.01528 g/mol
% molar mass of 1,1,1-Trichloroethane: 133.4 g/mol
% initial moles of water and NAPL
n_w = Fluid.rho_w*V_w / 18.01528;     % moles of water
n_n = Fluid.rho_n*V_n / 131.4;      % moles of NAPL

% volume of NAPL and water vapor (initially 0 b/c there's none)
V_gw_tot = zeros(size(T));
V_gn_tot = zeros(size(T));

n_gw_tot = zeros(size(T));
n_gn_tot = zeros(size(T));

times = [t];
temps = [mean(mean(T - 273.15*ones(Grid.Nz,Grid.Nx)))];

% w = VideoWriter('coboil_cells_unif_perm_Sn025.avi')
% open(w);

% v = VideoWriter('gasFlow_unif_perm.avi')
% open(v);
old_T = T;
% co_boil = zeros(size(T));

Tdata = [T(1,1)];

Tn = T;
Tn1 = temp_v3(Grid, T, Q, lambda, heat_cap, f_l, f_r);

% cb_time = [];
% cb_temp = [];
% 

% define left and right extractors (this is subject to change as cells
% dry out)
extractor_l = cell(Grid.Nz, 1);
extractor_r = cell(Grid.Nz, 1);
for i = 1:Grid.Nz
    extractor_l{i,1} = [i, 1];
    extractor_r{i,1} = [i, Grid.Nx];
end
    
Sg_vals = [];   t_Sg = [];

cb = 0;

co_boil = zeros(size(T));
old_cb = zeros(size(T));

cb_w_T = zeros(size(T));
cb_n_T = zeros(size(T));
cb_wn_T = zeros(size(T));

[A, rhsvec] = compute_A(Grid, lambda,...
    heat_cap, f_l, f_r);


% [L,U] = lu(A);

% i = 1;

%%

while t < (2*1296000)
    
%     if cb == 864000
%         Grid.dt = 360;
%         
%     elseif cb == 1036800
%         Grid.dt = 180;
%     end
    
    t = t + Grid.dt;
     
    % compute temp
    T = temp_v7(Grid, T, Q,lambda, heat_cap, f_l, f_r, co_boil);

% BDF2
%     if cb == 0
%         T_0 = reshape(Tn', Grid.Nx*Grid.Nz, 1);    
%         T_1 = reshape(Tn1', Grid.Nx*Grid.Nz, 1);
%         T = reshape(T', Grid.Nx*Grid.Nz, 1);
% 
%         T = A \ (2*T_1 - 0.5*T_0 - rhsvec);
%         
%         T = reshape(T, Grid.Nx, Grid.Nz)';
%     else
%         
%          T = temp_v6(Grid, Tn1, Tn, Q, lambda,...
%            heat_cap, f_l, f_r);
%        
%     end
% 
%     Tn = Tn1;    Tn1 = T;
    
    Tdata = [Tdata; T(1,1)];
    
    % compute vapor pressure using Antoine eqn
    % EDIT: P_nv = 0 if no NAPL is present in the cell
    P_nv = exp(19.796*ones(size(T)) - 2289.8*ones(size(T))...
        ./(T - 83.445*ones(size(T)))) .* (S_n ~= 0);
    P_wv = exp(23.195*ones(size(T)) - 3814*ones(size(T))...
        ./(T - 46.29*ones(size(T)))) .* (S_w > Fluid.S_r);
   
   
    co_boil = ((P_wv + P_nv) >= (P_w + Fluid.P_D));
    
%     co_boil = co_boil + ((P_wv + P_nv) >= (P_w + Fluid.P_D));
    
%     co_boil = co_boil > 0;
    
%     figure(4)
%     colormap([1 1 1; 0 0 1]);
%     image((co_boil) .* 255);
    
%     frame =  getframe(gcf);
%     writeVideo(w, frame);
        
    % check if T reaches co-boiling temp
    if any(co_boil, 'all') == 1
        
        cb = 1;
        
        
        
%         if co_boil(i,1) == 1
%             cb_time = [cb_time; t];
%             cb_temp = [cb_temp; T(i,1)];
%             i = i+1; z 0.5 dz0.05 x20 dx 0.05
%         end
%         
         % Compute Q which is given by Q = div(lambda * grad(T))
         
         % Gradient of T
         % since MATLAB's gradient solver uses a one sided
         % difference at the boundary it's not very accurate. But since
         % we know the heat flux at the boundary, we can edit the boundary
         % to be the heat flux defined earlier.
%          [T_x, T_z] = gradient(T, x, z);
%          T_x(:,1) = f_l;    T_x(:,end) = f_r;
%          T_z(1,:) = 0;      T_z(end,:) = 0;

         T_l = T(:,2) - 2*f_l*Grid.dx;      
         T_r = T(:,end-1) + 2*f_r*Grid.dx;
         
         T1 = [T_l, T, T_r];
         
         T_t = T1(2,:);  T_b = T1(end-1,:);
         
         T1 = [T_t; T1; T_b];
         
         [T_x, T_z] = gradient(T1, [-Grid.dx, x, Grid.dx + Grid.x],...
             [-Grid.dz, z, Grid.dz + Grid.z]);
         
         T_x = T_x(2:end-1, 2:end-1);       T_z = T_z(2:end-1, 2:end-1);
       
         Q = divergence(x, z, lambda.*T_x, lambda.*T_z);
                  
         % Q > 0 only in cells that reach co-boiling temperature, otherwise
         % Q = 0
         Q = Q .* co_boil;
%          Q = abs(Q);
         
%          Q = Q .* (Q > 0);
         Q_p = Q .* (Q > 0);
         
         % compute the moles of vapor produced using the energy
         % balance equation and Dalton's law Grid.dx*Grid.dz*Grid.dt
         n_gw = ((Q_p*Grid.dx*Grid.dz*Grid.dt)./ (Fluid.L_n*(P_wv./P_nv) +...
             Fluid.L_w*ones(size(Q))));                 % energy balance
         n_gn = n_gw.*P_nv./P_wv;                     % Dalton's law
         
         n_gn = n_gn .* co_boil;        n_gw = n_gw .* co_boil;
         n_gn(isnan(n_gn)) = 0;         n_gw(isnan(n_gw)) = 0;
         
         % Compute capillary and gas pressures
         % EDIT: there's a separate function for computing P_c and P_g
         P_c = Fluid.P_D .* ((max(S_w,Fluid.S_r) - Fluid.S_r*...
             ones(size(S_w)))./((1-Fluid.S_r)...
             *ones(size(S_w))-S_n)).^(-1/Fluid.Lambda);

         % local gas pressure   
         P_g =  P_c + P_w;

         % compute volume of gas using the ideal gas law
         V_gn = (8.314462*(n_gn.*T) ./ P_g);       % NAPL vapor
         
         if min(min(V_n - V_gn)) < 0
            V_gn =  V_gn .* (V_n > V_gn) + V_n .* (V_n < V_gn);
         end
         
         V_gn_tot = V_gn_tot + V_gn;
         n_gn_tot = n_gn_tot + V_gn.*P_g ./ (8.314462*T);
         
         V_gw = (8.314462*(n_gw.*T) ./ P_g);       % water vapor
         
         % S_w can only be in the range [S_r, 1] so the volume of water can
         % only be within the range for those S_w values
         if min(min(V_w - V_gw)) < 0
            V_gw = V_gw .* (V_w > V_gw) + V_w .* (V_w < V_gw);
         end
         
         V_gw_tot = V_gw_tot + V_gw;
         n_gw_tot = n_gw_tot + V_gw.*P_g ./ (8.314462*T);
         
         V_n = (V_n - V_gn) .* (V_n >= V_gn);
         V_w = (V_w - V_gw) .* (V_w >= V_gw);
      
         % update saturations
         S_n = V_n / V_cell;
         S_w = (V_w / V_cell) + Fluid.S_r;
%          S_g = (V_gw_tot + V_gn_tot)/V_cell;
         S_g = S_g + (V_gw + V_gn)/V_cell;
         
         % Compute capillary and gas pressures
         % capillary pressure
%          P_c = Fluid.P_D .* ((max(S_w,Fluid.S_r) - Fluid.S_r*...
%              ones(size(S_w)))./((1-Fluid.S_r)...
%              *ones(size(S_w))-S_n)).^(-1/Fluid.Lambda);
% 
%          % local gas pressure   
%          P_g =  P_c + P_w;        
         
%          % macro-IP
         if max(max(S_g)) > Fluid.S_gcr  
%              
             figure(3)
             colormap([1 1 1; 0 0 1]);
             image((S_g > Fluid.S_gcr) .* 255);
             
%              frame =  getframe(gcf);
%              writeVideo(v, frame);

             [S_g, S_w, S_n, Q, T] = macroIP(S_g, S_n, S_w, P_w, Q, T,...
                 co_boil, Fluid);
             
             figure(3)
             colormap([1 1 1; 0 0 1]);
             image((S_g > Fluid.S_gcr) .* 255);
             
             
%              
% %              frame =  getframe(gcf);
% %              writeVideo(v, frame);
%              
%              % we will need to recompute the volume of water since it will
%              % move around in the cell
%              V_w = (S_w - Fluid.S_r) * V_cell;
%              
%              Sg_vals = [Sg_vals; S_g(:,1)'];   t_Sg = [t_Sg; t];
             
         end
         
         % INCOMPLETE
         % Remove NAPL and water vapor from system
         % Since gas gets moved around in MIP it's hard to keep track of
         % the volume of NAPL and water vapor separately. But one way is to
         % keep track of how much liquid water and NAPL is still left in
         % the system. 
         
         % check for dried out regions adjacent to extractors
         for i = 1:Grid.Nz
           
             % left extractor
             if S_w(extractor_l{i,1}(1,1), extractor_l{i,1}(1,2) + 1)...
                     == Fluid.S_r
                 
                 extractor_l{i,1}(1,2) = extractor_l{i,1}(1,2) + 1;
                 
             end
             
             % right extractor
             if S_w(extractor_r{i,1}(1,1), extractor_r{i,1}(1,2) - 1)...
                     == Fluid.S_r
                 
                 extractor_r{i,1}(1,2) = extractor_r{i,1}(1,2) - 1;
                 
             end
             
         end
         
         % update moles of water and NAPL
         n_w = Fluid.rho_w*V_w / 18.01528;     % moles of water
         n_n = Fluid.rho_n*V_n / 131.4;      % moles of NAPL
              
         % update heat capacity
         heat_cap = S_w*Fluid.por*Fluid.rho_w*Fluid.C_pw + ...
             S_n*Fluid.por*Fluid.rho_n*Fluid.C_pn + ...
             (1-Fluid.por)*Fluid.rho_s*Fluid.C_ps;

         % update thermal conductivity
         K_e = (kappa*(1 - S_g))./(1 + (kappa - 1)*(1 - S_g));    
         lambda = K_e*(lambda_sat - lambda_dry) + lambda_dry;
         
         Q_old = Q;
    end

%     t = t + Grid.dt;
    
    old_T = T;
    
    times = [times; t];
    temps = [temps; mean(mean(T - 273.15*ones(Grid.Nz,Grid.Nx)))];
end

% close(v);
% close(w);

T = T - 273.15*ones(Grid.Nz,Grid.Nx);

% figure(1)
% [cs, hs] = contourf(xx,zz,flip(T,1));
% set(hs,'EdgeColor','none')
% colorbar
% caxis([10 600])
% colormap(jet)
% set(gca, 'Fontsize', 20)

% figure(2)
% plot(times, temps, 'Linewidth', 4)
% xlabel('t (seconds)')
% ylabel('Temperature (celsius)')
% set(gca, 'Fontsize', 20)
% title('Average temperature')
