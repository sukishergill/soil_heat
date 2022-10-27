function Q = compute_Q(T, f_l, f_r, lambda, Grid)

dx = Grid.dx; dz = Grid.dz;
Nx = Grid.Nx; Nz = Grid.Nz;

lambda = lambda';

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

           asup2(idx) = (l_jm+l_jp)/dz^2;
           
           [asub1(idx),asup1(idx), rhsvec(idx)] = ...
               solve_i(i, j, lambda, Nx, dx, f_l, f_r);
           
           
       elseif j == Nz
           
           l_jp = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i,j-1))/...
               (3*lambda(i,j)-lambda(i,j-1));
           l_jm = 2*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+lambda(i,j-1));

           asub2(idx) = ((l_jm+l_jp)/dz^2);
           
           [asub1(idx),asup1(idx),rhsvec(idx)] = ...
               solve_i(i, j, lambda, Nx, dx, f_l, f_r);
                                      
       else
           l_jp = 2*(lambda(i,j)*lambda(i,j+1))/(lambda(i,j)+lambda(i,j+1));
           l_jm = 2*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+lambda(i,j-1));  

           asub2(idx) = (l_jm/dz^2);
           asup2(idx) = (l_jp/dz^2);
           
           [asub1(idx),asup1(idx), rhsvec(idx)] = ...
               solve_i(i, j, lambda, Nx, dx, f_l, f_r);
            
       end
       
       % main diagonal 
       a1(idx) = -1*(asub1(idx)+asub2(idx)+asup1(idx)+asup2(idx));

       idx = idx + 1;
   end
end

asup1 = [zeros(1,1); asup1(1:end-1)];
asup2 = [zeros(Nx,1); asup2(1:end-Nx);];

asub1 = [asub1(2:end); zeros(1,1)];
asub2 = [asub2(Nx+1:end); zeros(Nx,1)];

A = spdiags( [ asub2, asub1, a1, asup1, asup2],...
    [-Nx, -1, 0, 1, Nx], Nx*Nz, Nx*Nz);


T = reshape(T', Nx*Nz, 1);

Q = A*T + rhsvec;

Q = reshape(Q, Nx, Nz)';

end

function [sub1, sup1, rhs] = solve_i(i1, j1, lambda,...
    Nx, dx, f_l, f_r)

if i1 == 1
   l_ip = 2*(lambda(i1,j1)*lambda(i1+1,j1))/(lambda(i1,j1)+lambda(i1+1,j1));
   l_im = (4*lambda(i1,j1)^2-2*lambda(i1,j1)*lambda(i1+1,j1))/...
       (3*lambda(i1,j1)-lambda(i1+1,j1));

   sup1 = (l_im+l_ip)/dx^2;
   sub1 = 0;
   rhs = -2*(l_im/dx^2)*dx*f_l;

elseif i1 == Nx
   l_ip = (4*lambda(i1,j1)^2-2*lambda(i1,j1)*lambda(i1-1,j1))/...
   (3*lambda(i1,j1)-lambda(i1-1,j1));
   l_im = 2*(lambda(i1,j1)*lambda(i1-1,j1))/(lambda(i1,j1)+lambda(i1-1,j1));

   sup1= 0;
   sub1= ((l_im+l_ip)/dx^2);
   rhs = 2*(l_ip/dx^2)*dx*f_r;
   
else
   l_ip = 2*(lambda(i1,j1)*lambda(i1+1,j1))/(lambda(i1,j1)+lambda(i1+1,j1));
   l_im = 2*(lambda(i1,j1)*lambda(i1-1,j1))/(lambda(i1,j1)+lambda(i1-1,j1));

   sub1 = (l_im/dx^2);
   sup1 = (l_ip/dx^2);
   rhs = 0;
end
end