function Tnew = temp_rad_1D(Grid, T, Q, lambda, ...
    heat_cap, f_R, f_inf)

dr = Grid.dr;   Nr = Grid.Nr;
dt = Grid.dt;

Q = dt * (Q ./ heat_cap);

% define diagonals
a1 = zeros(Nr, 1);          % main diagonal

asub = zeros(Nr, 1);        % subdiagonal
asup = zeros(Nr, 1);        % superdiagonal

rhsvec = zeros(Nr, 1);      % vector for the boundary terms

for j = 1:Nr
    
    r_jm = (j-0.15)*dr;
    r_jp = (j+0.5)*dr;
    
    if j == 1
        
        rl_m = r_jm*(4*lambda(j)^2 - 2*lambda(j)*lambda(j+1))...
            / (3*lambda(j) - lambda(j+1));
%         rl_m = (2*j*(j-1)*lambda(j)*(2*lambda(j) - lambda(j+1))) / ...
%             (j*lambda(j) + (j-1)*(2*lambda(j) - lambda(j+1)));
        
        rl_p = r_jp*(2*lambda(j)*lambda(j+1))/(lambda(j) + lambda(j+1));
%         rl_p = (2*j*(j+1)*lambda(j)*lambda(j+1)) / ...
%             (j*lambda(j) + (j+1)*lambda(j+1));
        
        a1(j) = 1 + (dt/(heat_cap(j)*dr^2*j))*(rl_p);
        asup(j) = -(dt*rl_p) / (heat_cap(j)*dr^2*j);
        
        rhsvec(j) = -dr*rl_m*f_R;
        
    elseif j == Nr
        
         rl_m = r_jm*(2*lambda(j)*lambda(j-1))/(lambda(j) + lambda(j-1));
%         rl_m = (2*j*(j-1)*lambda(j)*lambda(j-1)) / ...
%             (j*lambda(j) + (j-1)*lambda(j-1));
        
        rl_p = r_jp*(4*lambda(j)^2 - 2*lambda(j)*lambda(j-1))...
            / (3*lambda(j) - lambda(j-1));
%         rl_p = (2*j*(j+1)*lambda(j)*(2*lambda(j)-lambda(j-1))) / ...
%             (j*lambda(j) + (j+1)*(2*lambda(j) - lambda(j-1)));

        
        a1(j) = 1 + (dt/(heat_cap(j)*dr^2*j))*(rl_m);
        asub(j) = -(dt*rl_m) / (heat_cap(j)*dr^2*j);
        
        rhsvec(j) = dr*rl_p*f_inf;
        
    else
        
       rl_m = r_jm*(2*lambda(j)*lambda(j-1))/(lambda(j) + lambda(j-1));
%         rl_m = (2*j*(j-1)*lambda(j)*lambda(j-1)) / ...
%             (j*lambda(j) + (j-1)*lambda(j-1));

       rl_p = r_jp*(2*lambda(j)*lambda(j+1))/(lambda(j) + lambda(j+1));
%         rl_p = (2*j*(j+1)*lambda(j)*lambda(j+1)) / ...
%             (j*lambda(j) + (j+1)*lambda(j+1));      
        
        a1(j) = 1 + (dt/(heat_cap(j)*dr^2*j))*(rl_m + rl_p);
        asub(j) = -(dt*rl_m) / (heat_cap(j)*dr^2*j);
        asup(j) = -(dt*rl_p) / (heat_cap(j)*dr^2*j);
        
    end
    
end

asub = [asub(2:end); zeros(1,1)];
asup = [zeros(1,1); asup(1:end-1)];

A = spdiags([ asub, a1, asup], [-1, 0, 1], Nr, Nr);

Tnew = A \ (T + rhsvec - Q);

end