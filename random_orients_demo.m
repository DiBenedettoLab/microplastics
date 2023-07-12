% random orientation demo

% random points on unit sphere using equal area projection method
% https://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
th = 2*pi*rand(10000,1); % random angle between 0 and 2pi
z = 2*(rand(10000,1) - 0.5); % random elevation between -1 and 1

x = sqrt(1-z.^2).*cos(th);
y = sqrt(1-z.^2).*sin(th);

figure; scatter3(x,y,z,'k.'); daspect([1 1 1]);

figure; histogram(x,20)
figure; histogram(y,20)
figure; histogram(z,20)

%%
% random points
x = 2*(rand(10000,1) - 0.5);
y = 2*(rand(10000,1) - 0.5);
z = 2*(rand(10000,1) - 0.5);

figure; scatter3(x,y,z,'k.'); daspect([1 1 1]);

figure; histogram(x,20)
figure; histogram(y,20)
figure; histogram(z,20)
