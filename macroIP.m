function [S_g, S_w, S_n, n_gn, n_gw] = macroIP_v8(S_g, S_n, S_w,...
    P_w, T, n_gw, n_gn, co_boil, V_cell, Fluid)

S_gcr = Fluid.S_gcr;

MIP_cells = (S_g > S_gcr);

[clusters, lw, num] = findClusters(MIP_cells);

% compute the open volume in the domain (volume not occupied by water and
% NAPL)
V_open = V_cell*(1 - (S_w + S_n));

% V_open = V_cell * S_g;

% compute capillary pressure
P_c = compute_Pc(P_w, S_w, S_n, Fluid);
P_g = P_c + P_w;

T_avg = avgCluster(T, clusters, lw, num);


P_g = compute_Pg(P_g, n_gn + n_gw, 8.314462, T_avg, V_open, MIP_cells,...
    lw, num);

% Redistrubute NAPL and water vapor
[n_gn, n_gw] = redistribute(n_gn, n_gw, P_g, V_open, T_avg, lw, num);
% n_w = n_tot - (n_n + n_gn + n_gw);

% P_c = P_g - P_w;

% S_w = ((P_c ./ Fluid.P_D) ^ (-Fluid.Lambda)) .* (1 - (Fluid.S_r + S_n))...
%     + Fluid.S_r;


[T_e, T_t] = threshold(P_c, P_w, S_w);

old_MIP = MIP_cells;

MIP_end = 0;

while MIP_end == 0
    
    % expansion
    [MIP_cells, S_g, S_w] = expand(S_g, S_n, S_w, P_g, T_e, co_boil,...
        clusters, MIP_cells, S_gcr);
    [clusters, lw, num] = findClusters(MIP_cells);

    V_open = V_cell*(1 - (S_w + S_n));

    % compute capillary pressure
    P_c = compute_Pc(P_w, S_w, S_n, Fluid);
    P_g = P_c + P_w;

    T_avg = avgCluster(T, clusters, lw, num);

    P_g = compute_Pg(P_g, n_gn + n_gw, 8.314462, T_avg, V_open, MIP_cells,...
        lw, num);

    % Redistrubute NAPL and water vapor
    [n_gn, n_gw] = redistribute(n_gn, n_gw, P_g, V_open, T_avg, lw, num);

    [T_e, T_t] = threshold(P_c, P_w, S_w);

    
    % mobilization
    [MIP_cells, S_g, S_w] = mobilize(S_g, S_n, S_w, T, T_e,...
        T_t, P_g ,co_boil, clusters, MIP_cells, S_gcr, lw, num);
    [clusters, lw, num] = findClusters(MIP_cells);

    V_open = V_cell*(1 - S_w - S_n);

    V_open = V_cell*(1 - (S_w + S_n));

    % compute capillary pressure
    P_c = compute_Pc(P_w, S_w, S_n, Fluid);
    P_g = P_c + P_w;

    T_avg = avgCluster(T, clusters, lw, num);

    P_g = compute_Pg(P_g, n_gn + n_gw, 8.314462, T_avg, V_open, MIP_cells,...
        lw, num);

    % Redistrubute NAPL and water vapor
    [n_gn, n_gw] = redistribute(n_gn, n_gw, P_g, V_open, T_avg, lw, num);

    [T_e, T_t] = threshold(P_c, P_w, S_w);

    
    % check if expansion and mobilization/fragmentation can be repeated
    
    if any(MIP_cells ~= old_MIP, 'all') == 1
        MIP_end = 0;
        old_MIP = MIP_cells;
        
    else
        MIP_end = 1;
        
    end
    
    
end

end