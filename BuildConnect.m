function M = BuildConnect(nx,ny)

ex=ones(nx,1);
Lx=spdiags([ex ex ex],[-1 0 1],nx,nx);
ey=ones(ny,1);
Ly=spdiags([ey ey],[-1 1],ny,ny);
 
Ix=speye(nx);
Iy=speye(ny);
M=kron(Iy,Lx)+kron(Ly,Ix);