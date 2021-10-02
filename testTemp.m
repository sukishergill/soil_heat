clear all; clc;

Grid.x = 5;     Grid.dx = 0.2;      Grid.Nx = Grid.x/Grid.dx;
Grid.z = 5;     Grid.dz = 0.25;      Grid.Nz = Grid.z/Grid.dz;


Grid.dt = 0.1;

T = 10*ones(Grid.Nz, Grid.Nx);

K_e = (1.9*(1 - zeros(Grid.Nz,Grid.Nx)))./...
    (1 + (1.9-1)*(1-zeros(Grid.Nz,Grid.Nx)));

lambda = K_e*(2.75 - 0.15) + 0.15;

heat_cap = ones(Grid.Nz,Grid.Nx)*0.3*1*4.184 + ...
    ones(Grid.Nz,Grid.Nx)*0.3*1.46*0.958 + ...
    (1 - 0.3)*2.7*0.8;

Q = zeros(Grid.Nz,Grid.Nx);

f_ext = 0;
f_int = 0;
T = temp(Grid, T, Q, lambda, heat_cap, f_ext, f_int);

%Anew = full(A);