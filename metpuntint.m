function [xk,iter,fx,z,mu] = metpuntint(funobj, fung,x0)

% Aproxima una solucion del problema
% Min f(x)
% Sujeto a g(x) ? 0.
% donde las funciones,f: R^n --> R  , g:R^n --> R^p
% por medio del metodo de puntos interiores y
% y actualizaci?n cuasi-Newton a la matriz hessiana del lagrangiano. 
% Input:
%       funobj: cadena de caracteres con el nombre del codigo de Matlab
%               con la expresion algebraica de f(x).
%       fung:   cadena de caracteres con el nombre del codigo de Matlab con la expresion 
%               algebraica de g(x).
%       x0:     es un punto inicial para el proceso de optimizacion que
%               satisface, g(x0) > 0
% Output:
%       xk:     aproximacion al minimo del problema.
%       fx:     valor de la funcion objetivo en el optimo xk.
%       iter:   numero de iteraciones en el proceso. 
%-----------------------------------------------------
% Funciones auxiliares:
%       gradiente.m: Aproxima el gradiente
%       jacobiana1.m: Aproxima la jacobiana
%       fun_obj.m: Funcion Objetivo
%       fun_restricciones.m: Restricciones del problema
%-----------------------------------------------------

% ======================================
%   Definicion de Parametros iniciales
% ======================================
    xk = x0;
    x_min = x0;
    gx = feval(fung,x0);
    p = length(gx);
    n = length(x0);

% ======================================
%   Parametros iniciales y valores dentro de la funcion:
% ======================================
    maxiter = 200; % Limite de Iteraciones.
    iter = 0; % Contador de Iteraciones.
    tol = 1e-05; % Tolerancia para cortar las condiciones de Karush-Kuhn-Tucker.
    sigma = 0.2; % Parametro de actualizacion para la matriz Hessiana.
    mu = ones(p,1); % Vector de multiplicadores de Lagrange.
    z = ones(p,1); % Vector de holgura inicial igual al vector de unos.
    B = eye(n); % Matriz inicial para la Hessiana del Lagrangeano en x0.
    gradf = gradiente(funobj,x0); % Gradiente de la funci??n Objetivo en x0.
    A = jacobiana1(fung,x0); % Jacobiana de la funci??n Objetivo evaluada en x0.

% ======================================
%   Definicion de variables auxiliares en el m??todo.
% ======================================
    v1 = gradf - A'*mu;
    v2 = -gx + z;
    v3 = mu.*z;
    r = 0.2;
    Z = diag(z);
    U = diag(mu);
    Z_inv = Z\(eye(length(z)));
    cnpo_agg = [];
    fun_obj_agg = [];

% ======================================
%   Plot del punto inicial
% ======================================
    subplot(2,2,[1 2])
    title('Poligonos');
    pbaspect([1 1 1])
    axis([-0.5 1 0 1]);
    hold on
    axis equal
    
% ======================================
%   Condiciones de primer orden
% ======================================
    cnpo = [v1;v2;v3]'; 
    
% ======================================
%   Proceso de optimizacion 
% ======================================
    while(norm(cnpo) > tol  && iter < maxiter)
 
% ======================================
%   Guardamos la norma de las condiciones de primer orden
%   y el valor de la funcion objetivo
% ======================================
    cnpo_agg(iter + 1) = norm(cnpo);
    fun_obj_agg(iter + 1) = -1 * feval(funobj,xk);
% ======================================
%       Sistema reducido
% ======================================      
        M = B + A'*Z_inv*U*A;
        ld = v1 - A'*Z_inv*U*v2 + A'*Z_inv*(v3 - r*ones(length(gx),1));

% ======================================
%       Solucion para Delta(w)
% ======================================
        Dx = -M\ld; 
        Dz = -v2 + A*Dx;
        Dmu = Z_inv*(-U*Z*ones(length(gx),1) + r*ones(length(gx),1) - U*Dz);

% ======================================
%       Recorte del paso (alfa)
% ======================================
        bt = []; gm = [];
        for k =1:p
            if (Dmu(k) < 0)
                bt = [bt; -(mu(k)/Dmu(k))];
            end
            if(Dz(k) < 0)
                gm = [gm; -(z(k)/Dz(k))];
            end
        end
        
        if (sum(bt==0))
            alfaz=1; 
        else 
            alfaz=min(bt);
        end
        
        if (sum(gm==0))
            alfamu=1; 
        else 
            alfamu=min(gm);
        end
        alpha=min(alfaz,alfamu);
        alpha=0.995*min([1 alpha]);
        
% ======================================
%       Actualizacion del punto objetivo
% ======================================
        xk = xk + alpha*Dx;

% ======================================
%       Actualizaci??n de variables para la siguiente iteraci??n
% ======================================
        mu = mu + alpha*Dmu;
        z = z + alpha*Dz;
        gradf = gradiente(funobj,xk);
        gx = feval(fung,xk);
        A = jacobiana1(fung,xk);
        v1 = gradf - A'*mu;
        v2 = -gx + z;
        v3 = mu.*z;
        r = 0.9*r;
        cnpo = [v1;v2;v3]';
        Z = diag(z);
        U = diag(mu);
        Z_inv = Z\(eye(length(z)));

% ======================================
%       Elementos de actualizaci??n de BFGS
% ======================================
        s = alpha*Dx;
        y = (gradf-A'*mu) - (gradiente(funobj,xk-alpha*Dx) - (jacobiana1(fung,xk-alpha*Dx)'*mu));
        if ( (dot(s,y)) >= sigma*(dot(s,B*s)) )
            lambda = 1;
        else
            lambda = ((1-sigma)*(dot(s,B*s)))/(dot(s,B*s) - dot(s,y));
        end
% ======================================
%       Funcion r(lambda)
% ======================================
        rv = (lambda - 1)*y + lambda*B*s;       

% ======================================  
%       Actualizacion de la Matriz Hessiana
% ======================================        
        B = B - (((B*s)*(s'*B))/(dot(s,B*s))) + ((dot(rv,rv))/(dot(rv,s)));
        
        if (cond(B)>10^-5)
            B = eye(length(xk));
        end
        
% ======================================
%       Pasamos a la siguiente iteraci??n
% ======================================
        iter = iter + 1;   

% ======================================
%       Graficar punto inicial
% ======================================
        if(iter == 1)
            [u,v]=pol2cart(xk((n/2)+1:n), xk(1:n/2));
            plot([u;0;u(1)], [v;0;v(1)],'-*r');
        end    
    end
    
% ======================================
%   Graficamos el resultado del m??todo
% ======================================
    [u,v]=pol2cart(xk((n/2)+1:n), xk(1:n/2));
    plot([u;0;u(1)], [v;0;v(1)],'-*g'); % Graficar punto final
% ======================================
%   Graficamos el resultado de fmincon
% ======================================    
    n=length(xk)/2;
    m=length(xk);
    LB=zeros(2*(n),1);
    UB=[ones((n),1); ones((n),1)*pi]; 
    options=optimoptions('fmincon','Display','off');
    [x2, ~] = fmincon('fun_obj',x_min,[],[],[],[],LB,UB,'fung_mincon',options); 
    theta2=x2(n+1:m);
    r2=x2(1:n);
    [x,y]=pol2cart(theta2, r2);
    plot([x;0;x(1)], [y;0;y(1)],'-ob')
    legend({'Poligono Inicial', 'Poligono Final', 'Poligono de Fminicon'},'Location','northwest')
    fx = feval(funobj,xk);
    hold off
    subplot(2,2,3)
    plot(cnpo_agg, '-ob')
    title('Condiciones de Primer Orden');
    legend({'Norma'}, 'Location', 'northwest')
    subplot(2,2,4)
    plot(fun_obj_agg, '-ob')
    title('Evoluci??n de la funcion objetivo');
    legend({'Area'}, 'Location', 'northwest')

end