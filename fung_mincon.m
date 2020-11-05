function [gx,hx] = fung_mincon(x)

% ======================================  
% Funcion que define las restricciones g(x) y h(x)
% Adaptada al metodo de Matlab fmincon
% ======================================  
n=length(x)/2;
r=x(1:n);
theta=x((n+1):2*n);
g1=[];
for i = 1:(n-1)
    g1=[g1;r(1+i:n).^2+r(1:n-i).^2-2*r(1+i:n).*r(1:n-i).*cos(theta(1+i:n)-theta(1:n-i))-1];
end
g2=theta(1:n-1)-theta(2:n);
gx=[g1;g2];
hx=[];