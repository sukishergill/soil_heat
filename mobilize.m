function [MIP_cells, S_g, S_w ,P_g] =  mobilize(S_g, S_n, S_w, T, T_e,...
    T_t, P_g ,co_boil, clusters, MIP_cells, S_gcr, lw, num, extractors)

Nx = size(S_g, 2);          Nz = size(S_g, 1);
S_w = S_w;

% Average temp across cluster
T = avgCluster(T, clusters, lw, num);

for i = 1:num
    
    % find all boundary clusters that are adjacent to cluster
    clust_bound = findAdjacent(clusters{i,1}, Nx, Nz);
    
    for j = 1:size(clusters{i,1}, 1)
        
        % Need to check each cell adjacent to cluster
        
        for l = 1:size(clust_bound,1)  
            
            % need to check if adjacent cell is not in another cluster
            if (MIP_cells(clust_bound(l,1), clust_bound(l,2)) == 0) && ...
                    (co_boil(clust_bound(l,1), clust_bound(l,2)) == 1)

                % condition for mobilization/fragmentation
                if (T_t(clusters{i,1}(j,1), clusters{i,1}(j,2)) > ...
                        T_e(clust_bound(l,1), clust_bound(l,2))) && ...
                        (extractors(clust_bound(l,1), clust_bound(l,2)) == 0)

                    % move gas from imbibed cell to newly gas saturated cell
                    S_g(clust_bound(l,1), clust_bound(l,2)) = ...
                        S_g(clust_bound(l,1), clust_bound(l,2)) + ...
                        (S_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) - 0.14);
    
%                     % Adjust the gas pressure in invaded cell using the
%                     % ideal gas law
%                     P_g(clust_bound(l,1), clust_bound(l,2)) = ...
%                         (P_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) * ...
%                         S_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) * ...
%                         T(clust_bound(l,1), clust_bound(l,2)) / ...
%                         (T(clusters{i,1}(j,1), clusters{i,1}(j,2))) * S_gcr);

                    % gas saturation in invaded cell is increased to S_gcr
                    S_g(clust_bound(l,1), clust_bound(l,2)) = S_gcr;

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
                    MIP_cells(clusters{i,1}(j,1), clusters{i,1}(j,2)) = 0;

                    % "replace" the cell that used to be gas occupied
                    clusters{i,1}(j,:) = clust_bound(l,:);

                    % remove cell from cluster boundary
                    clust_bound(l,:) = [];

                    % check to see if clusters were fragmented or collided
                    % use an if statement to see if size of original cluster is
                    % the same as the new cluster.

                    % break out of inner most for loop because only one cell at
                    % a time can be imbibed
                    break

                end
            end
        end
    end
end

end