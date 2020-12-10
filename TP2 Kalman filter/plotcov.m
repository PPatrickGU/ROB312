function plotcov(x,P,style)
NP=50;
alpha  = 2*pi/NP*(0:NP);
circle = [cos(alpha);sin(alpha)];
ns = 1; 
C = chol(P(1:2,1:2))';
ellip = ns*C*circle;
X = x(1) + ellip(1,:);
Y = x(2) + ellip(2,:);
plot(X,Y,style)
