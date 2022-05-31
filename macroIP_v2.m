function [S_g, S_w, S_n] = macroIP_v2(S_g, S_n, S_w, P_g, P_c,...
    P_w, S_gcr, Grid)


[T_e, T_t] = threshold(P_c, P_w, S_w);

MIP_cells = (S_g > S_gcr);
non_MIP = (S_g <= S_gcr);

clusters = findCluster(MIP_cells);

old_clusters = clusters;

MIP_end = 0;

while MIP_end == 0 
    % equilibriate
    [S_g, S_w, P_g, P_c] =  equilibriate(S_g, S_w, P_g, P_c,... 
        MIP_cells, non_MIP);
    
    % expansion
    [clusters, MIP_cells, non_MIP] = expand(P_g, T_e, clusters,...
        MIP_cells, non_MIP, Grid);    
    [S_g, S_w, P_g, P_c] =  equilibriate(S_g, S_w, P_g, P_c,... 
        MIP_cells, non_MIP);
    [T_e, T_t] = threshold(P_c, P_w, S_w);
    
    % mobilization
    
    MIP_cells = (S_g > S_gcr);      non_MIP = (S_g <= S_gcr);
    
    % This is to ensure that after averaging, the average S_g in the 
    % cluster doesn't fall below the critical gas saturation
    if any(MIP_cells(:) == 1) ~= 1
        break
        
    else
        % We will need to recompute the clusters in case any clusters
        % split or combined 
        clusters = findClusters(MIP_cells);
    end
    
    % check if expansion and mobilization/fragmentation can be repeated
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
    
    old_clusters = clusters;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE THRESHOLDS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T_e, T_t] = threshold(P_c, P_w, S_w)

% parameter based on capillary pressure measurements
% P_t/P_e ratio
alpha = 0.57;

% entry and terminal pressures
P_e = P_c .* S_w;       P_t = alpha * P_e;

% drainage and imbibition thesholds
T_e = P_e + P_w;        T_t = P_t + P_w;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% FIND AND BUILD CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clusters = findCluster(MIP_cells)

[lw, num] = bwlabel(MIP_cells, 4);

clusters = cell(num, 1);

for i = 1:num
    [r,c] = find(lw == i);
    clusters{i,1} = [r,c];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% EQUILIBRIATE VALUES OVER CLUSTERS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S_g, S_w, P_g, P_c] =  equilibriate(S_g, S_w, P_g, P_c,...
    MIP_cells, non_MIP)

S_g = S_g .* non_MIP + mean2(nonzeros(S_g .* MIP_cells)) .* MIP_cells;
S_w = S_w .* non_MIP + mean2(nonzeros(S_w .* MIP_cells)) .* MIP_cells;
P_g = P_g .* non_MIP + mean2(nonzeros(P_g .* MIP_cells)) .* MIP_cells;
P_c = P_c .* non_MIP + mean2(nonzeros(P_c .* MIP_cells)) .* MIP_cells;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXPANSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clusters, MIP_cells, non_MIP] = expand(P_g, T_e, clusters,...
    MIP_cells, non_MIP, Grid)

for i = 1:size(clusters,1)
    
    for j = 1:size(clusters{i,1}, 1)
       % N cell
       % check the N cell isn't already in the cluster
       if ismember( [clusters{i,1}(j,1)-1, clusters{i,1}(j,2)],...
               clusters{i,1},'rows') == 0
           
          % check N cell isn't a boundary cell
          if (clusters{i,1}(j,1)-1) ~= 0
              
              % condition for expansion
              if P_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) > ...
                      T_e(clusters{i,1}(j,1)-1, clusters{i,1}(j,2))
                 
                  clusters{i,1} = [clusters{i,1} ;...
                      clusters{i,1}(j,1)-1, clusters{i,1}(j,2)];
                  
                  MIP_cells(clusters{i,1}(j,1)-1, clusters{i,1}(j,2)) = 1;
                  non_MIP(clusters{i,1}(j,1)-1, clusters{i,1}(j,2)) = 0;

              end
          end
          
       end
       
       % S cell
       % check the S cell isn't already in the cluster
       if ismember( [clusters{i,1}(1,1)+1, clusters{i,1}(j,2)],...
               clusters{i,1},'rows') == 0
           
          % check S cell isn't a boundary cell
          if (clusters{i,1}(1,1)+1) ~= Grid.Nz + 1
              
              % condition for expansion
              if P_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) > ...
                      T_e(clusters{i,1}(j,1)+1, clusters{i,1}(j,2))
                 
                  clusters{i,1} = [clusters{i,1} ;...
                      clusters{i,1}(j,1)+1, clusters{i,1}(j,2)];
                  
                  MIP_cells(clusters{i,1}(j,1)+1, clusters{i,1}(j,2)) = 1;
                  non_MIP(clusters{i,1}(j,1)+1, clusters{i,1}(j,2)) = 0;
                  
              end              
              
          end
          
       end
       
       % E cell
       % check E cell isn't already in the cluster
       if ismember( [clusters{i,1}(j,1), clusters{i,1}(j,2)+1],...
               clusters{i,1},'rows') == 0
           
          % check E cell isn't a boundary cell 
          if (clusters{i,1}(j,2)+1) ~= Grid.Nx + 1
              
              % condition for expansion
              if P_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) > ...
                      T_e(clusters{i,1}(j,1), clusters{i,1}(j,2)+1)
                 
                  clusters{i,1} = [clusters{i,1} ;...
                      clusters{i,1}(j,1), clusters{i,1}(j,2)+1];
                  
                  MIP_cells(clusters{i,1}(j,1), clusters{i,1}(j,2)+1) = 1;
                  non_MIP(clusters{i,1}(j,1), clusters{i,1}(j,2)+1) = 0;
  
              end               
          end
              
              
       end
       
       % W cell
       % check if W cell isn't already in the cluster
       if ismember( [clusters{i,1}(j,1), clusters{i,1}(j,2)-1],...
               clusters{i,1},'rows') == 0
           
          % check for expansion and/or mobilization
          if (clusters{i,1}(j,2)-1) ~= 0
              
              % condition for expansion
              if P_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) > ...
                      T_e(clusters{i,1}(j,1), clusters{i,1}(j,2)-1)
                 
                  clusters{i,1} = [clusters{i,1} ;...
                      clusters{i,1}(j,1), clusters{i,1}(j,2)-1];
                  
                  MIP_cells(clusters{i,1}(j,1), clusters{i,1}(j,2)-1) = 1;
                  non_MIP(clusters{i,1}(j,1), clusters{i,1}(j,2)-1) = 0;         
              end               
          end
       end
        
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MOBILIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INOMPLETE
function [clusters, MIP_cells, non_MIP] = mobilize(T_e, T_t, clusters,...
    boundary_clusters, MIP_cells, non_MIP, Grid)


for i = 1:size(clusters,1)
    
    for j = 1:size(clusters{i,1}, 1)
        
        % N cell
        % check if the N cell isn't already in the cluster
        if ismember( [clusters{i,1}(j,1)-1, clusters{i,1}(j,2)],...
               clusters{i,1},'rows') == 0
           
           % check if N cell isn't a boundary cell
           if (clusters{i,1}(j,1)-1) ~= 0
               
           end
        end
    end
end
end