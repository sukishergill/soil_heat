function [P_g, P_c] = computePressure(P_w, S_w, S_n, Fluid)

P_D = Fluid.P_D;    S_r = Fluid.S_r;    Lambda = Fluid.Lambda;

P_c = P_D * ((max(S_w, S_r) - S_r*ones(size(S_w)))...
    ./((1 - S_r)*ones(size(S_w)) - S_n)).^(-1/Lambda);

P_g = P_c + P_w;

end