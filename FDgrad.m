% function [fx, fy] = FDgrad(F,Fxa,Fxb,dx)
%
% Calculate the 2D gradient of a function on a uniform grid.
%
% The function is assumed ordered as the result of meshgrid.
%
% The boundaries are: fx(a,y) = fxa, fx(b,y) = fxb
% fy(x,c) = 0, fy(x,d) = 0
%
% The grid spacing is passed as dx where dx is a scalar or
% a vector of [dx; dy]
%
%
function [fx, fy] = FDgrad(F, Fxa, Fxb, dx)

if length(dx) == 2
    dy2 = 0.5/dx(2);
    dx2 = 0.5/dx(1);
    
else
    dx2 = 0.5/dx;
    dy2 = dx2;
end

[M,N] = size(F);
fx = 0*F;
fy = 0*F;

I = 2:(M-1);
J = 2:(N-1);

fx(:,J) = (F(:, J+1) - F(:, J-1))*dx2;
fy(I,:) = (F(I+1, :) - F(I-1, :))*dy2;

fx(:,1) = Fxa;
fx(:,N) = Fxb;

end