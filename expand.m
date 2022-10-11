function [MIP_cells, S_g, S_w] = expand(S_g, S_n, S_w, P_g, T_e,...
    co_boil, clusters, MIP_cells, S_gcr)

Nx = size(S_g, 2);          Nz = size(S_g, 1);

for i = 1:size(clusters,1)
    
    for j = 1:size(clusters{i,1}, 1)
       % N cell
       % check the N cell isn't already in the cluster
       if ismember( [clusters{i,1}(j,1)-1, clusters{i,1}(j,2)],...
               clusters{i,1},'rows') == 0
           
          % check N cell isn't a boundary cell
          if ((clusters{i,1}(j,1)-1) ~= 0) && ... 
                  (co_boil(clusters{i,1}(j,1)-1, clusters{i,1}(j,2)) == 1)
              
              % condition for expansion
              if P_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) > ...
                      T_e(clusters{i,1}(j,1)-1, clusters{i,1}(j,2))
                 
                  clusters{i,1} = [clusters{i,1} ;...
                      clusters{i,1}(j,1)-1, clusters{i,1}(j,2)];
                  
                  S_g(clusters{i,1}(j,1)-1, clusters{i,1}(j,2)) = S_gcr;
                  S_w(clusters{i,1}(j,1)-1, clusters{i,1}(j,2)) = 1 - ...
                      (S_g(clusters{i,1}(j,1)-1, clusters{i,1}(j,2)) + ...
                      S_n(clusters{i,1}(j,1)-1, clusters{i,1}(j,2)));
                  MIP_cells(clusters{i,1}(j,1)-1, clusters{i,1}(j,2)) = 1;
%                   non_MIP(clusters{i,1}(j,1)-1, clusters{i,1}(j,2)) = 0;

              end
          end
          
       end
       
       % S cell
       % check the S cell isn't already in the cluster
       if ismember( [clusters{i,1}(j,1)+1, clusters{i,1}(j,2)],...
               clusters{i,1},'rows') == 0
           
          % check S cell isn't a boundary cell
          if ((clusters{i,1}(j,1)+1) ~= Nz + 1) && ... 
                  (co_boil(clusters{i,1}(j,1)+1, clusters{i,1}(j,2)) == 1)
              
              % condition for expansion
              if P_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) > ...
                      T_e(clusters{i,1}(j,1)+1, clusters{i,1}(j,2))
                 
                  clusters{i,1} = [clusters{i,1} ;...
                      clusters{i,1}(j,1)+1, clusters{i,1}(j,2)];
                  
                  S_g(clusters{i,1}(j,1)+1, clusters{i,1}(j,2)) = 0.15;
                  MIP_cells(clusters{i,1}(j,1)+1, clusters{i,1}(j,2)) = 1;
%                   non_MIP(clusters{i,1}(j,1)+1, clusters{i,1}(j,2)) = 0;
                  
              end              
              
          end
          
       end
       
       % E cell
       % check E cell isn't already in the cluster
       if ismember( [clusters{i,1}(j,1), clusters{i,1}(j,2)+1],...
               clusters{i,1},'rows') == 0
           
          % check E cell isn't a boundary cell 
          if ((clusters{i,1}(j,2)+1) ~= Nx + 1) && ... 
                  (co_boil(clusters{i,1}(j,1), clusters{i,1}(j,2)+1) == 1)
              
              % condition for expansion
              if P_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) > ...
                      T_e(clusters{i,1}(j,1), clusters{i,1}(j,2)+1)
                 
                  clusters{i,1} = [clusters{i,1} ;...
                      clusters{i,1}(j,1), clusters{i,1}(j,2)+1];
                  
                  S_g(clusters{i,1}(j,1), clusters{i,1}(j,2)+1) = 0.15;
                  MIP_cells(clusters{i,1}(j,1), clusters{i,1}(j,2)+1) = 1;
%                   non_MIP(clusters{i,1}(j,1), clusters{i,1}(j,2)+1) = 0;
  
              end               
          end
              
              
       end
       
       % W cell
       % check if W cell isn't already in the cluster
       if ismember( [clusters{i,1}(j,1), clusters{i,1}(j,2)-1],...
               clusters{i,1},'rows') == 0
           
          % check for expansion and/or mobilization
          if ((clusters{i,1}(j,2)-1) ~= 0) && ... 
                  (co_boil(clusters{i,1}(j,1), clusters{i,1}(j,2)-1) == 1)
              
              % condition for expansion
              if P_g(clusters{i,1}(j,1), clusters{i,1}(j,2)) > ...
                      T_e(clusters{i,1}(j,1), clusters{i,1}(j,2)-1)
                 
                  clusters{i,1} = [clusters{i,1} ;...
                      clusters{i,1}(j,1), clusters{i,1}(j,2)-1];
                  
                  S_g(clusters{i,1}(j,1), clusters{i,1}(j,2)-1) = 0.15;
                  MIP_cells(clusters{i,1}(j,1), clusters{i,1}(j,2)-1) = 1;
%                   non_MIP(clusters{i,1}(j,1), clusters{i,1}(j,2)-1) = 0;         
              end               
          end
       end
        
    end
end


end
