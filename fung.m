function [gx] = fung(x)


% ======================================  
% Funcion que define las restricciones g(x)
% ======================================  
n=length(x)/2;
r=x(1:n);
theta=x((n+1):2*n);
g1=zeros(n*(n-1),1);
k=1;
for i = 1:(n-1)
    for j = (i+1):n
        g1(k)=1-r(i)^2-r(j)^2+2*r(i)*r(j)*cos(theta(i)-theta(j));
        k=k+1;
    end
end

g2=theta(2:n)-theta(1:n-1);
g3=theta;
g4=pi-theta;
g5=r;
g6=1-r;

gx=[g1;g2;g3;g4;g5;g6];