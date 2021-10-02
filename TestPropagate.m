nx = 200;
ny = 300;
[x,y] = meshgrid((0:(nx-1))/(nx-1),(0:(ny-1))/(ny-1));
U = abs(sin(3*x-y).*cos(16*(y-1/2).^2));

thresh = .5;

subplot(2,1,1);
contourf(x,y,U,[.1 .2 .5]);

Clusters = propagate(U, thresh);

subplot(2,1,2);
contourf(x,y,U,[.1 .2 .5]);
hold on;
N = size(Clusters);
for i=1:N(2)
    plot(x(Clusters(:,i)==1),y(Clusters(:,i)==1),'.','MarkerSize',20);
end