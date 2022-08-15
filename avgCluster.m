function B =  avgCluster(A, clusters, lw, num_clust)

B = zeros(size(A));

for i = 1:num_clust
    
    B = B + avg_val(A, clusters{i,1}) .* (lw == i);
    
end

B = B + A .* (lw == 0);

end

function m = avg_val(M, cluster)

sum_clust = 0;

for j = 1:size(cluster, 1)
    sum_clust = sum_clust + M(cluster(j,1), cluster(j,2));
end

m = sum_clust/size(cluster, 1);

end