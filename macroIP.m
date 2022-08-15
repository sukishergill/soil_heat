function [S_g, S_w, S_n, Q, T] = macroIP(S_g, S_n, S_w, P_w, Q, T,...
    co_boil, Fluid, Grid)

S_gcr = Fluid.S_gcr;

MIP_cells = (S_g > S_gcr);
% non_MIP = (S_g <= S_gcr);     % might not need this

[clusters, lw, num] = findClusters(MIP_cells);

% [S_g, S_w, S_n] = avgCluster(S_g, S_w, S_n, clusters, lw, num);
% [P_g, P_c] = computePressure(P_w, S_w, S_n, Fluid);
% [T_e, T_t] = threshold(P_c, P_w, S_w);

S_g = avgCluster(S_g, clusters, lw, num);
S_w = 1 - (S_g + S_n);

[P_g, P_c] = computePressure(P_w, S_w, S_n, Fluid);
P_g = avgCluster(P_g, clusters, lw, num);

[T_e, T_t] = threshold(P_c, P_w, S_w);

% old_clusters = clusters;
old_MIP = MIP_cells;

MIP_end = 0;

while MIP_end == 0
    
    % expansion
    [MIP_cells, S_g, S_w] = expand(S_g, S_n, S_w, P_g, T_e, co_boil,...
        clusters, MIP_cells, Grid, S_gcr);
    [clusters, lw, num] = findClusters(MIP_cells);

    S_g = avgCluster(S_g, clusters, lw, num);
    S_w = 1 - (S_g + S_n);

    [P_g, P_c] = computePressure(P_w, S_w, S_n, Fluid);
    P_g = avgCluster(P_g, clusters, lw, num);

    [T_e, T_t] = threshold(P_c, P_w, S_w);
    
    % mobilization
    [MIP_cells, S_g, S_w] = mobilize(S_g, S_n, S_w, T_e, T_t, co_boil,...
        clusters, MIP_cells, Grid, S_gcr);
    [clusters, lw, num] = findClusters(MIP_cells);

    S_g = avgCluster(S_g, clusters, lw, num);
    S_w = 1 - (S_g + S_n);
    
    [P_g, P_c] = computePressure(P_w, S_w, S_n, Fluid);
    P_g = avgCluster(P_g, clusters, lw, num);
    
    [T_e, T_t] = threshold(P_c, P_w, S_w);
    
    % check if expansion and mobilization/fragmentation can be repeated
    
    if any(MIP_cells ~= old_MIP, 'all') == 1
        MIP_end = 0;
        old_MIP = MIP_cells;
        
    else
        MIP_end = 1;
        
    end
    
%     % Check to see if the number of clusters has changed which will
%     % indicate if any clusters fragmented or collided
%     if size(clusters, 1) ~= size(old_clusters, 1)
%         MIP_end = 0;
%         old_clusters = clusters;
%         
%     % If the cluster size has not changed then we look at if the size of 
%     % each cluster has changed. Cluster sizes will change if it has 
%     % undergone expansion
%     
%     else
%         for i = 1:size(clusters, 1)
%         
%             if size(clusters{i,1}, 1) ~= size(old_clusters{i,1}, 1)
%                 MIP_end = 0;
%                 old_clusters = clusters;
%                 break
%                 
%             else
%                 
%                 MIP_end = 1;
%                 
%                 for j = 1:size(clusters{i,1}, 1)
%                     
%                     if ismember([clusters{i,1}(j,1), clusters{i,1}(j,2)],...
%                             old_clusters{i,1}, 'rows') == 0
%                         
%                         MIP_end = 1;
%                         break
%                     
%                     end            
% 
%                 end
%             end
%         end
%     end
    
end

% Q = avgCluster(Q, clusters, lw, num);
% T = avgCluster(T, clusters, lw, num);

end