function P_g = compute_Pg(P_g, n_g, R, T, V_open, MIP_cells, lw, num)

% n_g = n_gn + n_gw;

n_g_tot = zeros(size(n_g));
V_tot = zeros(size(V_open));

for i = 1:num
    
    n_g_tot = n_g_tot + sum(sum(n_g .* (lw == i))) .* (lw == i);
    V_tot = V_tot + sum(sum(V_open .* (lw == i))) .* (lw == i);
    
end

P_avg = (R*n_g_tot.*T) ./ V_tot;
P_avg(isnan(P_avg)) = 0;

P_g = P_avg .* MIP_cells + P_g .* (MIP_cells == 0);

end