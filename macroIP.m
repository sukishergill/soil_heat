function Sg = macroIP(S_w, S_g, S_n, Fluid, Grid, T)

% parameter based on capillary pressure measurements
% P_t/P_e ratio
alpha = 0.57;               

% define pressures
P_e = Fluid.P_D*S_w;                         % entry pressure
P_t = alpha*P_e;                             % terminal pressure
P_w = Fluid.rho_w*9.8*Grid.z + 1.01*10^5;    % water pressure

% gas pressure
P_g = ((S_w - Fluid.S_gcr)./(1 - Fluid.S_gcr - S_n)).^...
    ((-1/Fluid.Lambda)*ones(size(T)))*Fluid.P_D + P_w*ones(size(T));

% compute thresholds
T_e = P_e + P_w;            % drainage threshold
T_t = P_t + P_w;            % imbibition threshold

% compute moles of NAPL vapor
[n_gn, n_gw] = NAPL_mole(Fluid, T);

% convert moles of vapor to volume using the ideal gas law (PV = nRT)
% !!!!!!!!!!!!!!!!!!!
% INCOMPLETE
% !!!!!!!!!!!!!!!!!!! 
V_n = (n_gn.*T) ./ P_g;         % NAPL vapor
V_w = (n_gw.*T) ./ P_g;         % water vapor

V = V_n + V_w;

% !!!!!!!!!!!!!!!!!!!
% INCOMPLETE
% !!!!!!!!!!!!!!!!!!!
% Convert V to S_g by figuring out the volume of the block and how much 
% space is available space.

%Sg = propagate(V, Fluid.S_gcr);

Sg = V/0.0015;


% logical array for saturation values greater that the critical gas
% saturation
%S_g2 = S_g >= Fluid.S_gcr;

end