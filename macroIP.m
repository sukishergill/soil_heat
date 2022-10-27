function [S_g, S_w, S_n, n_gn, n_gw] = macroIP(S_g, S_n, S_w,...
    P_w, Q, T, V_gw, V_gn, n_gw, n_gn, co_boil, V_cell, Fluid)


S_gcr = Fluid.S_gcr;

MIP_cells = (S_g > S_gcr);

[clusters, lw, num] = findClusters(MIP_cells);

% compute the "open" volume in the domain
V_open = V_cell*(1 - S_w - S_n);

P_g = (8.314462*(n_gn + n_gw).*T) ./ V_open;
P_g = avgCluster(P_g, clusters, lw, num);

% Redistrubute NAPL and water vapor
[n_gn, n_gw] = redistribute(n_gn, P_g, V_open, T, lw, num);

P_c = computePressure(P_w, S_w, S_n, Fluid);

[T_e, T_t] = threshold(P_c, P_w, S_w);

old_MIP = MIP_cells;

MIP_end = 0;

while MIP_end == 0
    
    % expansion
    [MIP_cells, S_g, S_w] = expand(S_g, S_n, S_w, P_g, T_e, co_boil,...
        clusters, MIP_cells, S_gcr);
    [clusters, lw, num] = findClusters(MIP_cells);

    V_open = V_cell*(1 - S_w - S_n);

    P_g = (8.314462*(n_gn + n_gw).*T) ./ V_open;
    P_g = avgCluster(P_g, clusters, lw, num);

    % Redistrubute NAPL and water vapor
    [n_gn, n_gw] = redistribute(n_gn, P_g, V_open, T, lw, num);

    P_c = computePressure(P_w, S_w, S_n, Fluid);

    [T_e, T_t] = threshold(P_c, P_w, S_w);

    
    % mobilization
    [MIP_cells, S_g, S_w] = mobilize(S_g, S_n, S_w, T, T_e,...
        T_t, P_g ,co_boil, clusters, MIP_cells, S_gcr, lw, num);
    [clusters, lw, num] = findClusters(MIP_cells);

    V_open = V_cell*(1 - S_w - S_n);

    P_g = (8.314462*(n_gn + n_gw).*T) ./ V_open;
    P_g = avgCluster(P_g, clusters, lw, num);

    % Redistrubute NAPL and water vapor
    [n_gn, n_gw] = redistribute(n_gn, P_g, V_open, T, lw, num);

    P_c = computePressure(P_w, S_w, S_n, Fluid);

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