function Tnew = temp(Grid, T, Q, lambda,...
    heat_cap, f_ext, f_int)

f = Q;

lambda = lambda';
heat_cap = heat_cap';

dx = Grid.dx; dz = Grid.dz;
Nx = Grid.Nx; Nz = Grid.Nz;
dt = Grid.dt;

a1 = repmat(ones(Nx,1), Nz, 1);       % main diagonal

% subdiagonals
asub1 = repmat(ones(Nx,1), Nz, 1);
asub2 = repmat(ones(Nx,1), Nz, 1);

% superdiagonals
asup1 = repmat(ones(Nx,1), Nz, 1);
asup2 = repmat(ones(Nx,1), Nz, 1);

idx = 1;

rhsvec = zeros(Nx*Nz,1);

%test_vals = zeros(10,16);

for j = 1:Nz
   for i = 1:Nx
       
       % Thermal conductivities (harmonic averaging
       % use if statement to deal with the boundary
       % first-order (linear) extrapolation to define lambdas outside the
       % domain
       
       if i == 1 && j == 1
           % lambda_{i+1/2,j}
           l_ip = 2*(lambda(i,j)*lambda(i+1,j))/(lambda(i,j)+lambda(i+1,j));
           %lambda_{i-1/2,j}
           l_im = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i+1,j))/...
               (3*lambda(i,j)-lambda(i+1,j));
           %lambda_{i,j+1/2}
           l_jp = 2*(lambda(i,j)*lambda(i,j+1))/(lambda(i,j)+lambda(i,j+1));
           %lambda_{i,j-1/2}
           l_jm = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i,j+1))/...
               (3*lambda(i,j)-lambda(i,j+1));

           sub1 = (-dt/heat_cap(i,j)) * (l_im/dx^2);
           sub2 = (-dt/heat_cap(i,j)) * (l_jm/dz^2);              
           sup1 = (-dt/heat_cap(i,j)) * (l_ip/dx^2);
           sup2 = (-dt/heat_cap(i,j)) * (l_jp/dz^2);
           
           % subdiagonals
           asub1(idx) = asub1(idx) * 0;              
           asub2(idx) = 0;

           % superdiagonals
           asup1(idx) = (sub1 + sup1);               
           asup2(idx) = (sub2 + sup2);

           rhsvec(idx) = 2*(-sub1*dx*f(idx) - sup2*dz*f(idx));
       
       elseif i == 1 && j == Nz
           
           l_ip = 2*(lambda(i,j)*lambda(i+1,j))/(lambda(i,j)+lambda(i+1,j));
           %lambda_{i-1/2,j}
           l_im = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i+1,j))/...
               (3*lambda(i,j)-lambda(i+1,j));
           %lambda_{i,j+1/2}
           l_jp = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i,j-1))/...
               (3*lambda(i,j)-lambda(i,j-1));
           %lambda_{i,j-1/2}
           l_jm = 2*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+lambda(i,j-1));

           sub1 = (-dt/heat_cap(i,j)) * (l_im/dx^2);
           sub2 = (-dt/heat_cap(i,j)) * (l_jm/dz^2);               
           sup1 = (-dt/heat_cap(i,j)) * (l_ip/dx^2);
           sup2 = (-dt/heat_cap(i,j)) * (l_jp/dz^2);

           % subdiagonals
           asub1(idx) = 0;              
           asub2(idx) = (sup2 + sub2);

           % superdiagonals
           asup1(idx) = (sub1 + sup1);               
           asup2(idx) = 0;

           rhsvec(idx) = 2*(-sub1*dx*f(idx) + sub2*dz*f(idx));
           
       elseif i == 1 && (j > 1 && j < Nz)
           % lambda_{i+1/2,j}
           l_ip = 2*(lambda(i,j)*lambda(i+1,j))/(lambda(i,j)+lambda(i+1,j));
           %lambda_{i-1/2,j}
           l_im = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i+1,j))/...
               (3*lambda(i,j)-lambda(i+1,j));
           %lambda_{i,j+1/2}
           l_jp = 2*(lambda(i,j)*lambda(i,j+1))/(lambda(i,j)+lambda(i,j+1));
           %lambda_{i,j-1/2}
           l_jm = 2*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+lambda(i,j-1));

           sub1 = (-dt/heat_cap(i,j)) * (l_im/dx^2);
           sub2 = (-dt/heat_cap(i,j)) * (l_jm/dz^2);             
           sup1 = (-dt/heat_cap(i,j)) * (l_ip/dx^2);
           sup2 = (-dt/heat_cap(i,j)) * (l_jp/dz^2);

           % subdiagonals
           asub1(idx) = 0;              
           asub2(idx) = sub2;

           % superdiagonals
           asup1(idx) = (sub1 + sup1);               
           asup2(idx) = sup2;

           rhsvec(idx) = 2*(-sub1*dx*f(idx));
          
       elseif i == Nx && j == 1
           % lambda_{i+1/2,j}
           l_ip = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i-1,j))/...
               (3*lambda(i,j)-lambda(i-1,j));
           %lambda_{i-1/2,j}
           l_im = 2*(lambda(i,j)*lambda(i-1,j))/(lambda(i,j)+lambda(i-1,j));
           %lambda_{i,j+1/2}
           l_jp = 2*(lambda(i,j)*lambda(i,j+1))/(lambda(i,j)+lambda(i,j+1));
           %lambda_{i,j-1/2}
           l_jm = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i,j+1))/...
               (3*lambda(i,j)-lambda(i,j+1));

           sub1 = (-dt/heat_cap(i,j)) * (l_im/dx^2);
           sub2 = (-dt/heat_cap(i,j)) * (l_jm/dz^2);               
           sup1 = (-dt/heat_cap(i,j)) * (l_ip/dx^2);
           sup2 = (-dt/heat_cap(i,j)) * (l_jp/dz^2);

           asub1(idx) = (sub1 + sup1);              
           asub2(idx) = 0;

           % superdiagonals
           asup1(idx) = 0;               
           asup2(idx) = sub2 + sup2;

           rhsvec(idx) = 2*(sup1*dx*f(idx) - sup2*dx*f(idx));
           
       elseif i == Nx && j == Nz
           % lambda_{i+1/2,j}
           l_ip = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i-1,j))/...
               (3*lambda(i,j)-lambda(i-1,j));
           %lambda_{i-1/2,j}
           l_im = 2*(lambda(i,j)*lambda(i-1,j))/(lambda(i,j)+lambda(i-1,j));
           %lambda_{i,j+1/2}
           l_jp = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i,j-1))/...
               (3*lambda(i,j)-lambda(i,j-1));
           %lambda_{i,j-1/2}
           l_jm = 2*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+lambda(i,j-1));

           sub1 = (-dt/heat_cap(i,j)) * (l_im/dx^2);
           sub2 = (-dt/heat_cap(i,j)) * (l_jm/dz^2);             
           sup1 = (-dt/heat_cap(i,j)) * (l_ip/dx^2);
           sup2 = (-dt/heat_cap(i,j)) * (l_jp/dz^2);

           % subdiagonals
           asub1(idx) = (sub1 + sup1);              
           asub2(idx) = sub2 + sup2;

           % superdiagonals
           asup1(idx) = 0;               
           asup2(idx) = 0;

           rhsvec(idx) = 2*(sup1*dx*f(idx) + sub2*dx*f(idx));
           
       elseif i == Nx && (j > 1 && j < Nz)
           % lambda_{i+1/2,j}
           l_ip = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i-1,j))/...
               (3*lambda(i,j)-lambda(i-1,j));
           %lambda_{i-1/2,j}
           l_im = 2*(lambda(i,j)*lambda(i-1,j))/(lambda(i,j)+lambda(i-1,j));
           %lambda_{i,j+1/2}
           l_jp = 2*(lambda(i,j)*lambda(i,j+1))/(lambda(i,j)+lambda(i,j+1));
           %lambda_{i,j-1/2}
           l_jm = 2*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+lambda(i,j-1));

           sub1 = (-dt/heat_cap(i,j)) * (l_im/dx^2);
           sub2 = (-dt/heat_cap(i,j)) * (l_jm/dz^2);
           sup1 = (-dt/heat_cap(i,j)) * (l_ip/dx^2);
           sup2 = (-dt/heat_cap(i,j)) * (l_jp/dz^2);

           % subdiagonals
           asub1(idx) = (sub1 + sup1);              
           asub2(idx) = sub2;

           % superdiagonals
           asup1(idx) = 0;               
           asup2(idx) = sup2;

           rhsvec(idx) = 2*(sup1*dx*f(idx));
           
       elseif j == 1 && (i > 1 && i < Nx)
           % lambda_{i+1/2,j}
           l_ip = 2*(lambda(i,j)*lambda(i+1,j))/(lambda(i,j)+lambda(i+1,j));
           %lambda_{i-1/2,j}
           l_im = 2*(lambda(i,j)*lambda(i-1,j))/(lambda(i,j)+lambda(i-1,j));
           %lambda_{i,j+1/2}
           l_jp = 2*(lambda(i,j)*lambda(i,j+1))/(lambda(i,j)+lambda(i,j+1));
           %lambda_{i,j-1/2}
           l_jm = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i,j+1))/...
               (3*lambda(i,j)-lambda(i,j+1));

           sub1 = (-dt/heat_cap(i,j)) * (l_im/dx^2);
           sub2 = (-dt/heat_cap(i,j)) * (l_jm/dz^2);              
           sup1 = (-dt/heat_cap(i,j)) * (l_ip/dx^2);
           sup2 = (-dt/heat_cap(i,j)) * (l_jp/dz^2);

           % subdiagonals
           asub1(idx) = sub1;              
           asub2(idx) = 0;

           % superdiagonals
           asup1(idx) = sup1;               
           asup2(idx) = sub2 + sup2;

           rhsvec(idx) = 2*(-sup2*dz*f(idx));
           
       elseif j == Nz && (i ~= 1 || i ~= Nx)
           % lambda_{i+1/2,j}
           l_ip = 2*(lambda(i,j)*lambda(i+1,j))/(lambda(i,j)+lambda(i+1,j));
           %lambda_{i-1/2,j}
           l_im = 2*(lambda(i,j)*lambda(i-1,j))/(lambda(i,j)+lambda(i-1,j));
           %lambda_{i,j+1/2}
           l_jp = (4*lambda(i,j)^2-2*lambda(i,j)*lambda(i,j-1))/...
               (3*lambda(i,j)-lambda(i,j-1));
           %lambda_{i,j-1/2}
           l_jm = 2*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+lambda(i,j-1));

           sub1 = (-dt/heat_cap(i,j)) * (l_im/dx^2);
           sub2 = (-dt/heat_cap(i,j)) * (l_jm/dz^2);
           sup1 = (-dt/heat_cap(i,j)) * (l_ip/dx^2);
           sup2 = (-dt/heat_cap(i,j)) * (l_jp/dz^2);

           % subdiagonals
           asub1(idx) = sub1;              
           asub2(idx) = sub2 + sup2;

           % superdiagonals
           asup1(idx) = sup1;               
           asup2(idx) = 0;

           rhsvec(idx) = 2*(sub2*dz*f(idx));
           
       elseif (i ~= 1 && i ~= Nx) && (j ~= 1 && j ~= Nz)
           % lambda_{i+1/2,j}
           l_ip = 2*(lambda(i,j)*lambda(i+1,j))/(lambda(i,j)+lambda(i+1,j));
           %lambda_{i-1/2,j}
           l_im = 2*(lambda(i,j)*lambda(i-1,j))/(lambda(i,j)+lambda(i-1,j));
           %lambda_{i,j+1/2}
           l_jp = 2*(lambda(i,j)*lambda(i,j+1))/(lambda(i,j)+lambda(i,j+1));
           %lambda_{i,j-1/2}
           l_jm = 2*(lambda(i,j)*lambda(i,j-1))/(lambda(i,j)+lambda(i,j-1));  
           
           sub1 = (-dt/heat_cap(i,j)) * (l_im/dx^2);
           sub2 = (-dt/heat_cap(i,j)) * (l_jm/dz^2);
           sup1 = (-dt/heat_cap(i,j)) * (l_ip/dx^2);
           sup2 = (-dt/heat_cap(i,j)) * (l_jp/dz^2);      
           
           % subdiagonals
           asub1(idx) = sub1;      
           asub2(idx) = sub2;
           
           % superdiagonals
           asup1(idx) = sup1;           
           asup2(idx) = sup2;
       end           
       
       % main diagonal 
       a1(idx) = 1 + (dt/heat_cap(i,j))*((l_ip+l_im)/dx^2 + (l_jp+l_jm)/dz^2);
       
       %test_vals(:,idx) = [i j l_ip l_im l_jp l_jm asup1(idx)...
       %    asup2(idx) asub1(idx) asub2(idx)];
       
       idx = idx + 1;
   end
end

asup1 = [zeros(1,1); asup1(1:end-1)];
asup2 = [zeros(Nx,1); asup2(1:end-Nx);];

asub1 = [asub1(2:end); zeros(1,1)];
asub2 = [asub2(Nx+1:end); zeros(Nx,1)];

A = spdiags( [ asub2, asub1, a1, asup1, asup2],...
    [-Nx, -1, 0, 1, Nx], Nx*Nz, Nx*Nz);

%spy(A)

% A = spdiags( [a1, a1, a1],...
%     [-1, 0, 1], Nx*Nz, Nx*Nz);

Q = reshape(Q', Nx*Nz, 1);

T = reshape(T', Nx*Nz, 1);

Tnew = A \ (T - rhsvec - Q);

Tnew = reshape(Tnew, Nx, Nz)';


end