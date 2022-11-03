function [new_n_gn, new_n_gw] = redistribute(n_gn, n_gw, P_g, V, T,...
    lw, num)
    
new_n_gn = zeros(size(n_gn));

for i = 1:num

    A = (lw == i);

    n = concentration(A, n_gn, V);

    new_n_gn = new_n_gn + n * V .* A;

end

new_n_gn = new_n_gn + (lw == 0) .* n_gn;

% P_gw = P_g - (8.314462*new_n_gn.*T)./V;

% new_n_gw = (P_gw.*V) ./ (8.314462*T);

new_n_gw = ((P_g .* V) ./ (8.314462*T)) - new_n_gn;

new_n_gw = new_n_gw.*(lw > 0) + n_gw.*(lw == 0);

end

function n = concentration(B, n_gn, V_gn)

n = sum(sum(n_gn .* B)) / sum(sum(V_gn .* B));

end