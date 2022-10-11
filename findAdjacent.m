function adjacent_cells = findAdjacent(cluster, Nx, Nz)

A = sparse(cluster(:,1), cluster(:,2), ones(size(cluster,1), 1),...
    Nz, Nx);
    
M = zeros(Nz + 2, Nx + 2);

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