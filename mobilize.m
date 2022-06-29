function [MIP_cells, S_g, S_w] =  mobilize(S_g, S_n, S_w, T_e,...
    T_t, clusters, MIP_cells, Grid, S_gcr)

% IMCOMPLETE: Still need to recompute the cluster after a cell is
% invaded incase it leads to two cells colliding, or after a cell is
% imbibed clusters could fragment.

for i = 1:size(clusters,1)
    
    % find all boundary clusters that are adjacent to cluster
    clust_bound = findAdjacent(clusters{i,1}, Grid);
    
    for j = 1:size(clusters{i,1}, 1)
        
        % Need to check each cell adjacent to cluster
        
        for l = 1:size(clust_bound,1)  

            % condition for mobilization/fragmentation
            if T_t(clusters{i,1}(j,1), clusters{i,1}(j,2)) > ...
                    T_e(clust_bound(l,1), clust_bound(l,2))
                
%                 % move gas from imbibed cell to newly gas saturated cell
%                 S_g(clust_bound(l,1), clust_bound(l,2)) = ...
%                     S_g(clust_bound(l,1), clust_bound(l,2)) + ...
%                     (S_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) - 0.14);
                
                % gas saturation in invaded cell is increased to S_gcr
                S_g(clust_bound(l,1), clust_bound(l,2)) = 0.15;
                
                % Adjust S_g in imbibed cell to 0.14 to represent 
                % trapped gas
                S_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) = 0.14;
                
                % Adjust the water saturation so that the equilibrium
                % law is still satisfied
                S_w(clust_bound(l,1), clust_bound(l,2)) = ...
                    1 - (S_g(clust_bound(l,1), clust_bound(l,2))...
                    + S_n(clust_bound(l,1), clust_bound(l,2)));
                
                S_w(clusters{i,1}(j,1), clusters{i,1}(j,2)) = 1 - ...
                    (S_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) + ...
                    S_n(clusters{i,1}(j,1), clusters{i,1}(j,2)));
                
                MIP_cells(clust_bound(l,1), clust_bound(l,2)) = 1;
%                 non_MIP(clust_bound(l,1), clust_bound(l,2)) = 0;
                
                MIP_cells(clusters{i,1}(j,1), clusters{i,1}(j,2)) = 0;
%                 non_MIP(clusters{i,1}(j,1), clusters{i,1}(j,2)) = 1;
                
                % "replace" the cell that used to be gas occupied
                clusters{i,1}(j,:) = clust_bound(l,:);
                
                % remove cell from cluster boundary
                clust_bound(l,:) = [];
                
                % check to see if clusters were fragmented or collided
                % use an if statement to see if size of original cluster is
                % the same as the new cluster.
                
                % break out of inner most for loop because only one cell at
                % a time can me imbibed
                break

            end
        end
    end
end

end