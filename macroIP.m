function [S_g, S_w, S_n, Q, T] = macroIP(S_g, S_n, S_w, P_w, Q, T,...
    co_boil, Fluid)


S_gcr = Fluid.S_gcr;

MIP_cells = (S_g > S_gcr);

[clusters, lw, num] = findClusters(MIP_cells);


S_g = avgCluster(S_g, clusters, lw, num);
S_w = 1 - (S_g + S_n);

P_c = computePressure(P_w, S_w, S_n, Fluid);

P_g = P_c + P_w;
P_g = avgCluster(P_g, clusters, lw, num);

[T_e, T_t] = threshold(P_c, P_w, S_w);

old_MIP = MIP_cells;

MIP_end = 0;

while MIP_end == 0
    
    % expansion
    [MIP_cells, S_g, S_w] = expand(S_g, S_n, S_w, P_g, T_e, co_boil,...
        clusters, MIP_cells, S_gcr);
    [clusters, lw, num] = findClusters(MIP_cells);

    S_g = avgCluster(S_g, clusters, lw, num);
    S_w = 1 - (S_g + S_n);

    P_c = computePressure(P_w, S_w, S_n, Fluid);
    P_g = avgCluster(P_g, clusters, lw, num);

    [T_e, T_t] = threshold(P_c, P_w, S_w);
    
    % mobilization
    [MIP_cells, S_g, S_w] = mobilize(S_g, S_n, S_w, T, T_e,...
        T_t, P_g ,co_boil, clusters, MIP_cells, S_gcr, lw, num);
    [clusters, lw, num] = findClusters(MIP_cells);

    S_g = avgCluster(S_g, clusters, lw, num);
    S_w = 1 - (S_g + S_n);
    
    P_c = computePressure(P_w, S_w, S_n, Fluid);
    P_g = avgCluster(P_g, clusters, lw, num);
    
    [T_e, T_t] = threshold(P_c, P_w, S_w);
    
    % check if expansion and mobilization/fragmentation can be repeated
    
    if any(MIP_cells ~= old_MIP, 'all') == 1
        MIP_end = 0;
        old_MIP = MIP_cells;
        
    else
        MIP_end = 1;
        
    end
    
    
end

% Q = avgCluster(Q, clusters, lw, num);
% T = avgCluster(T, clusters, lw, num);

end