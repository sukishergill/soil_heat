function [clusters, lw, num] = findClusters(MIP_cells)

[lw, num] = bwlabel(MIP_cells, 4);

clusters = cell(num, 1);

for i = 1:num
    [r,c] = find(lw == i);
    clusters{i,1} = [r,c];
end
end