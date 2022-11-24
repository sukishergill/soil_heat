function Tnew = temp_v6(Grid, T1, T0, Q, lambda,...
    heat_cap, f_l, f_r)


dx = Grid.dx; dz = Grid.dz;
Nx = Grid.Nx; Nz = Grid.Nz;
dt = Grid.dt;

Q = dt * (Q ./ heat_cap);
Q = reshape(Q', Nx*Nz, 1);

lambda = lambda';
heat_cap = heat_cap';
heat_cap = 2*heat_cap;

a1 = repmat(zeros(Nx,1), Nz, 1);       % main diagonal

% subdiagonals
asub1 = repmat(zeros(Nx,1), Nz, 1);
asub2 = repmat(zeros(Nx,1), Nz, 1);

% superdiagonals
asup1 = repmat(zeros(Nx,1), Nz, 1);
asup2 = repmat(zeros(Nx,1), Nz, 1);

idx = 1;

rhsvec = zeros(Nx*Nz,1);


for j = 1:Nz
   for i = 1:Nx
       
       if j == 1
           l_jp = 2*(lambda(i,j)*lambda(i,j+1))/(lambda(i,j)+lambda(i,j+1));
           l_jm = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i,j+1))/...
               (3*lambda(i,j)-lambda(i,j+1));

           asup2(idx) = (-dt/heat_cap(i,j)) * (l_jm+l_jp)/dz^2;
           
           [asub1(idx),asup1(idx), rhsvec(idx)] = ...
               solve_i(i, j, heat_cap, lambda, Nx, dx, dt, f_l(j), f_r(j));
           
           
       elseif j == Nz
           
           l_jp = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i,j-1))/...
               (3*lambda(i,j)-lambda(i,j-1));
           l_jm = 2*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+lambda(i,j-1));

           asub2(idx) = (-dt/heat_cap(i,j))*((l_jm+l_jp)/dz^2);
           
           [asub1(idx),asup1(idx),rhsvec(idx)] = ...
               solve_i(i, j, heat_cap, lambda, Nx, dx, dt, f_l(j), f_r(j));
                                      
       else
           l_jp = 2*(lambda(i,j)*lambda(i,j+1))/(lambda(i,j)+lambda(i,j+1));
           l_jm = 2*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+lambda(i,j-1));  

           asub2(idx) = (-dt/heat_cap(i,j)) * (l_jm/dz^2);
           asup2(idx) = (-dt/heat_cap(i,j)) * (l_jp/dz^2);
           
           [asub1(idx),asup1(idx), rhsvec(idx)] = ...
               solve_i(i, j, heat_cap, lambda, Nx, dx, dt, f_l(j), f_r(j));
            
       end
       
       % main diagonal 
       a1(idx) = 1.5 - (asub1(idx)+asub2(idx)+asup1(idx)+asup2(idx));

       idx = idx + 1;
   end
end

asup1 = [zeros(1,1); asup1(1:end-1)];
asup2 = [zeros(Nx,1); asup2(1:end-Nx);];

asub1 = [asub1(2:end); zeros(1,1)];
asub2 = [asub2(Nx+1:end); zeros(Nx,1)];

A = spdiags( [ asub2, asub1, a1, asup1, asup2],...
    [-Nx, -1, 0, 1, Nx], Nx*Nz, Nx*Nz);

A = (2/3)*A;

T0 = reshape(T0', Nx*Nz, 1);    T1 = reshape(T1', Nx*Nz, 1);

% [L, U] = lu(A);

% Tnew = U \ (L \ (2*T1 - 0.5*T0 - rhsvec - Q));

Tnew = A \ ((4/3)*T1 - (1/3)*T0 - rhsvec - Q);

Tnew = reshape(Tnew, Nx, Nz)';

end

function [sub1, sup1, rhs] = solve_i(i1, j1, heat_cap, lambda,...
    Nx, dx, dt, f_l, f_r)

if i1 == 1
   l_ip = 2*(lambda(i1,j1)*lambda(i1+1,j1))/(lambda(i1,j1)+lambda(i1+1,j1));
   l_im = (4*lambda(i1,j1)^2-2*lambda(i1,j1)*lambda(i1+1,j1))/...
       (3*lambda(i1,j1)-lambda(i1+1,j1));

   sup1 = (-dt/heat_cap(i1,j1)) * (l_im+l_ip)/dx^2;
   sub1 = 0;
   rhs = -2*(-dt/heat_cap(i1,j1))*(l_im/dx^2)*dx*f_l;

elseif i1 == Nx
   l_ip = (4*lambda(i1,j1)^2-2*lambda(i1,j1)*lambda(i1-1,j1))/...
   (3*lambda(i1,j1)-lambda(i1-1,j1));
   l_im = 2*(lambda(i1,j1)*lambda(i1-1,j1))/(lambda(i1,j1)+lambda(i1-1,j1));

   sup1= 0;
   sub1= (-dt/heat_cap(i1,j1))*((l_im+l_ip)/dx^2);
   rhs = 2*(-dt/heat_cap(i1,j1))*(l_ip/dx^2)*dx*f_r;
   
else
   l_ip = 2*(lambda(i1,j1)*lambda(i1+1,j1))/(lambda(i1,j1)+lambda(i1+1,j1));
   l_im = 2*(lambda(i1,j1)*lambda(i1-1,j1))/(lambda(i1,j1)+lambda(i1-1,j1));

   sub1 = (-dt/heat_cap(i1,j1)) * (l_im/dx^2);
   sup1 = (-dt/heat_cap(i1,j1)) * (l_ip/dx^2);
   rhs = 0;
end
end