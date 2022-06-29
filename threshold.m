function [T_e, T_t] = threshold(P_c, P_w, S_w)

% parameter based on capillary pressure measurements
% P_t/P_e ratio
alpha = 0.57;

% entry and terminal pressures
P_e = P_c .* S_w;       P_t = alpha * P_e;

% drainage and imbibition thesholds
T_e = P_e + P_w;        T_t = P_t + P_w;

end