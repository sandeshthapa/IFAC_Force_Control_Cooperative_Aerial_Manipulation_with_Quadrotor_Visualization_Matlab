function h = circle(x,y,r,n)
% Plots a circle
hold on
th = 0:n:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off