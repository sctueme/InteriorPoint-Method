function [fx] = fun_obj(x0)



% ======================================
% Funcion objetivo adaptada para minimizacion
% ====================================== 
n = length(x0)/2;
m = length(x0);
r=x0(1:n);
t=x0((n+1):m);
flag=r(2:n).*r(1:n-1).*sin(t(2:n)-t(1:n-1));
fx = sum(flag);
fx=-fx/2;
end