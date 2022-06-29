function [S_g, S_w, S_n] = macroIP(S_g, S_n, S_w, P_w, Fluid, Grid)

S_gcr = Fluid.S_gcr;

MIP_cells = (S_g > S_gcr);
% non_MIP = (S_g <= S_gcr);     % might not need this

[clusters, lw, num] = findClusters(MIP_cells);

[S_g, S_w, S_n] = avgCluster(S_g, S_w, S_n, clusters, lw, num);
[P_g, P_c] = computePressure(P_w, S_w, S_n, Fluid);
[T_e, T_t] = threshold(P_c, P_w, S_w);

old_clusters = clusters;

MIP_end = 0;

while MIP_end == 0
    
    % expansion
    [MIP_cells, S_g, S_w] = expand(S_g, S_n, S_w, P_g, T_e, clusters,...
        MIP_cells, Grid, S_gcr);
    [clusters, lw, num] = findClusters(MIP_cells);
    [S_g, S_w, S_n] = avgCluster(S_g, S_w, S_n, clusters, lw, num);
    [P_g, P_c] = computePressure(P_w, S_w, S_n, Fluid);
    [T_e, T_t] = threshold(P_c, P_w, S_w);
    
    % mobilization
%     [MIP_cells, S_g, S_w] = mobilize(S_g, S_n, S_w, T_e, T_t, clsuters,...
%         MIP_cells, Grid, S_gcr);
%     [clusters, lw, num] = findClusters(MIP_cells);
%     [S_g, S_w, S_n] = avgCluster(S_g, S_w, S_n, clusters, lw, num);
%     [P_g, P_c] = computePressure(P_w, S_w, S_n, Fluid);
%     [T_e, T_t] = threshold(P_c, P_w, S_w);
    
    % check if expansion and mobilization/fragmentation can be repeated
    % Needs to be edited b/c this looks at whether the cluster size changes
    % but in a mobilization event it doesn't change
    if size(clusters, 1) ~= size(old_clusters, 1)
        MIP_end = 0;
        old_clusters = clusters;
        
    else
        for i = 1:size(clusters, 1)
        
            if size(clusters{i,1}, 1) ~= size(old_clusters{i,1}, 1)
                MIP_end = 0;
                old_clusters = clusters;
                break
                
            else
                MIP_end = 1;
                
            end
        end
    end
    
end

end