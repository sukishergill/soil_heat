function [Sg, Sw, Sn] =  avgCluster(S_g, S_w, S_n,...
    clusters, lw, num_clust)

Sg = zeros(size(S_g));
Sw = Sg;    

% Still need to figure out if we need to average S_n over a cluster
% Sn = Sg;
Sn = S_n;

for i = 1:num_clust
    
    Sg = Sg + avg_val(S_g, clusters{i,1}, lw, i);
    Sw = Sw + avg_val(S_w, clusters{i,1}, lw, i);
%     Sn = Sn + avg_val(S_n, clusters{i,1}, lw, i);
    
end

Sg = S_g .* (lw == 0) + Sg;
Sw = S_w .* (lw == 0) + Sw;
% Sn = S_n .* (lw == 0) + Sn;

end

function M = avg_val(A, cluster, N, c)

sum_clust = 0;

for j = 1:size(cluster, 1)
    sum_clust = sum_clust + A(cluster(j,1), cluster(j,2));
end

M = (sum_clust/size(cluster, 1)) .* (N == c);

end