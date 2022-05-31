function Tnew = temp_rad_2D(Grid, T, Q, lambda, ...
    heat_cap, f_i, f_o)

dr = Grid.dr;       Nr = Grid.Nr;
dth = Grid.dth;     Nth = Grid.Nth;
dt = Grid.dt;

Q = dt * (Q ./ heat_cap);
Q = reshape(Q', Nr*Nth, 1);

lambda = lambda';
heat_cap = heat_cap';


% define diagonals
a1 = repmat(zeros(Nr,1), Nth, 1);       % main diagonal

% subdiagonals
asub1 = repmat(zeros(Nr,1), Nth, 1);
asub2 = repmat(zeros(Nr,1), Nth, 1); 

% superdiagonals
asup1 = repmat(zeros(Nr,1), Nth, 1);
asup2 = repmat(zeros(Nr,1), Nth, 1);

idx = 1;

rhsvec = zeros(Nr*Nth, 1);

for j = 1:Nth
   for i = 1:Nr
       
       if j == 1
           l_jp = 2*dr*i*(lambda(i,j) * lambda(i,j+1))/(lambda(i,j)+...
               lambda(i,j+1));
           l_jm = dr*i*(4*lambda(i,j)^2-2*lambda(i,j)*lambda(i,j+1))/...
               (3*lambda(i,j)-lambda(i,j+1));
           
           asup2(idx) = (-dt/heat_cap(i,j)) * (l_jm+l_jp)/(dth*i*dr)^2;
          
           
       elseif j == Nth
           
           l_jp = dr*i*(4*lambda(i,j)^2-2*lambda(i,j)*lambda(i,j-1))/...
               (3*lambda(i,j)-lambda(i,j-1));
           l_jm = 2*dr*i*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+...
               lambda(i,j-1));
           
           asub2(idx) = (-dt/heat_cap(i,j)) * (l_jm+l_jp)/(dth*i*dr)^2;
           
       else
           l_jp = 2*dr*i*(lambda(i,j) * lambda(i,j+1))/(lambda(i,j)+...
               lambda(i,j+1));          
           l_jm = 2*dr*i*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+...
               lambda(i,j-1));
           
           asub2(idx) = (-dt/heat_cap(i,j)) * (l_jm)/(dth*i*dr)^2;
           asup2(idx) = (-dt/heat_cap(i,j)) * (l_jp)/(dth*i*dr)^2;
       end
       
       [asub1(idx), asup1(idx), rhsvec(idx)] = ...
           solve_i(i, j, heat_cap, lambda, Nr, dr, dt, f_i, f_o);
       
       a1(idx) = 1 - (asub1(idx)+asub2(idx)+asup1(idx)+asup2(idx));
       
       idx = idx + 1;
       
   end
end

asup1 = [zeros(1,1); asup1(1:end-1)];
asup2 = [zeros(Nr,1); asup2(1:end-Nr);];

asub1 = [asub1(2:end); zeros(1,1)];
asub2 = [asub2(Nr+1:end); zeros(Nr,1)];

A = spdiags( [ asub2, asub1, a1, asup1, asup2],...
    [-Nr, -1, 0, 1, Nr], Nr*Nth, Nr*Nth);


T = reshape(T', Nr*Nth, 1);

Tnew = A \ (T - rhsvec - Q);

Tnew = reshape(Tnew, Nr, Nth)';


end

function [sub1, sup1, rhs] = solve_i(i1, j1, heat_cap, lambda,...
    Nr, dr, dt, f_i, f_o)

r_im = (i1 - 0.5)*dr;       r_ip = (i1 + 0.5)*dr;

if i1 == 1
    
    l_im = r_im*(4*lambda(i1,j1)^2 - 2*lambda(i1,j1)*lambda(i1+1,j1))...
        / (3*lambda(i1,j1) - lambda(i1+1,j1));   
    l_ip = r_ip*(2*lambda(i1,j1)*lambda(i1+1,j1))/(lambda(i1,j1)...
        + lambda(i1+1,j1));
    
    sup1 = (-dt/heat_cap(i1,j1)) * (l_im+l_ip)/(dr^3*i1);
    sub1 = 0;
    rhs = -2*(-dt/heat_cap(i1,j1))*(l_im/(dr^3*i1))*dr*f_i;
    
elseif i1 == Nr
    
    l_im = r_im*2*(lambda(i1,j1)*lambda(i1-1,j1))/(lambda(i1,j1)+...
        lambda(i1-1,j1));
    l_ip = r_ip*(4*lambda(i1,j1)^2-2*lambda(i1,j1)*lambda(i1-1,j1))/...
       (3*lambda(i1,j1)-lambda(i1-1,j1));
   
    sub1 = (-dt/heat_cap(i1,j1)) * (l_im+l_ip)/(dr^3*i);
    sup1 = 0;
    rhs = -2*(-dt/heat_cap(i1,j1))*(l_ip/(dr^3*i))*dr*f_o;
    
else
    l_im = r_im*2*(lambda(i1,j1)*lambda(i1-1,j1))/(lambda(i1,j1)+...
        lambda(i1-1,j1));
    l_ip = r_ip*(2*lambda(i1,j1)*lambda(i1+1,j1))/(lambda(i1,j1)...
        + lambda(i1+1,j1));
    
    sub1 = (-dt/heat_cap(i1,j1)) * (l_im)/(dr^3*i);
    sup1 = (-dt/heat_cap(i1,j1)) * (l_ip)/(dr^3*i1);
    rhs = 0;
end

end
