%Poligono de 9 vertices


np = 9;

nv=1:np-1;

%Particion del semicirculo

theta = (pi/np) * nv;
r = 0.2; %valor del radio
%Punto inicial
X0=r*[cos(theta);sin(theta)];


[theta0,r0]=cart2pol(X0(1,:),X0(2,:));

x0=[r0';theta0'];

funobj = 'fun_obj';
fung = 'fung';

[x, iter, fx,z,mu] = metpuntint(funobj,fung,x0);
areaf = -fx;

fprintf('El area final del poligono de %i vertices es %f\n', np,areaf)