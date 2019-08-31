clear all, close all, clc
X = [2 4 5 7];
Y = [0 2 1 0];
Z = [1 3 5 7];
r = 0.5
Coords = [X;Y;Z]
Coords_inc = Coords;
Coords_inc(:,3) = Coords(:,4) - Coords(:,3);
Coords_inc(:,2) = Coords(:,3) - Coords(:,2);
Coords_inc(:,1) = Coords(:,2) - Coords(:,1)
x = 1;
y = 2;
z = 3;
hold on
color = ['r' 'g' 'b' 'k'];
theta = -atan(Coords_inc(y,:)./Coords_inc(z,:)); % angle around x
phi = atan(Coords_inc(x,:)./Coords_inc(z,:));% angle 2 (y)
psi = atan(Coords_inc(y,:)./Coords_inc(x,:));% angle 3 (z)
for i=1:4
      Rx = [1 0 0; 0 cos(theta(i)) -sin(theta(i)); 0 sin(theta(i)) cos(theta(i))]
      Ry = [cos(phi(i)) 0 sin(phi(i)); 0 1 0; -sin(phi(i)) 0 cos(phi(i))]
      Rz = [cos(psi(i)) -sin(psi(i)) 0; sin(psi(i)) cos(psi(i)) 0; 0 0 1]
      t = 0:0.1:2*pi
      x = cos(t)
      y = sin(t)
      z = zeros(1,length(t))
      tol = [X(i); Y(i); Z(i)]
      base = [x;y;z]
      rotated = Rx*Ry*Rz*base
      plot3(rotated(1,:)+X(i),rotated(2,:)+Y(i),rotated(3,:)+Z(i), color(i))
      axis equal
end