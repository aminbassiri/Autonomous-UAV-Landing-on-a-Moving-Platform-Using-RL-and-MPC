function c = plotCircle(x0,y0,r)
theta = linspace(0,2*pi,100)';
c(:,1) = r*cos(theta) + x0;
c(:,2) = r*sin(theta) + y0;
end