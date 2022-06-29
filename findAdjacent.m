function adjacent_cells = findAdjacent(cluster, Grid)

A = sparse(cluster(:,1), cluster(:,2), ones(size(cluster,1), 1),...
    Grid.Nz, Grid.Nx);
    
M = zeros(Grid.Nz + 2, Grid.Nx + 2);

% shift left
M(2:end-1, 1:end-2) = M(2:end-1, 1:end-2) + A;

% shift right
M(2:end-1, 3:end) = M(2:end-1, 3:end) + A;

% shift up
M(1:end-2, 2:end-1) = M(1:end-2, 2:end-1) + A;

% shift down
M(3:end, 2:end-1) = M(3:end, 2:end-1) + A;

M = M(2:end-1, 2:end-1);

M = (M ~= 0);
M = M - A;

[r, c] = find(M == 1);

adjacent_cells = [r, c];
    
end