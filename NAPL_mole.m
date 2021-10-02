function [n_gn, n_gw] = NAPL_mole(Fluid, T)

% compute vapor pressure using Antoine eqn
P_wv = 10.^(8.07131*ones(size(T)) - 1730.63*ones(size(T))...
    ./(233.426*ones(size(T)) + T));
P_nv = 10.^(7.11*ones(size(T)) - 1367.05*ones(size(T))...
    ./(235.809*ones(size(T))+T));


% !!!!!!!!!!!!!!!!!!!
% INCOMPLETE
% !!!!!!!!!!!!!!!!!!!
% include integral term from equation

% moles of NAPL vapor
n_gn = ones(size(T)) ./ (Fluid.L_w*(P_wv./P_nv) + Fluid.L_n);

% moles of water vapor
n_gw = n_gn .* (P_wv./P_nv);

end