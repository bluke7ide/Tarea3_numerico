function sigmas = calculo_s2(xn, yn)
    %
    %   Calcula los sigmas para poder realizar s2
    %
    %   Inputs
    %       xn: coordenadas x de los puntos a interpolar
    %       yn: coordenadas y de los puntos a interpolar
    %
    %   Outputs
    %       sigmas: coeficientes de la función s2 

    m = length(xn);
    x = zeros(m);
    y = zeros(m,1);
    for i = 1:(m-1)
        % diferencias de xn
        x(i,i) = (xn(i+1)-xn(i))/2;
        x(i,i+1) = x(i,i);
        % diferencias de yn
        y(i) = yn(i+1) - yn(i);
    end
    % condición de sigmas, que el último debe de ser igual al primero
    x(end, 1) = 1;
    x(end, end) = -1;
    y(end) = 0;

    sigmas = x\y;
end

m = 15;
xn = sort(rand(m+1,1));
f = @(x) cos(2.*pi.*x);
xn(1) = 0;
xn(end) = 1;
yn = f(xn);
sigmas = calculo_s2(xn, yn);

function s2 = s2(sigmas, x, xn, yn)
    %
    %   Calcula la función s2 en un punto x
    %
    %   Inputs 
    %       sigmas: calculados con calculo_s2
    %       x: punto a evaluar
    %       xn: coordenadas x de los puntos
    %       yn: coordenadas y de los puntos
    %
    %   Outputs
    %       s2: valor de s2 en x 

    % busca en qué intervalo está
    nodes = find(xn>x);
    if x == 1 
        lower = length(xn)-1;
    else 
        lower = nodes(1)-1;
    end
    
    % y procede a hacer la función
    s1 = sigmas(lower);
    s2 = sigmas(lower+1);
    y1 = yn(lower);
    x1 = xn(lower);
    x2 = xn(lower + 1);
    s2 = y1 + s1.*(x - x1) + ((s2 - s1)/(2 * (x2 - x1)))*(x-x1)^2;
end
figure;
xx = linspace(0,1);
yy = arrayfun(@(x) s2(sigmas, x, xn, yn), xx);
plot(xx, f(xx),'', 'DisplayName', 'f(x)')
hold on;
plot(xx, yy, 'r--', 'DisplayName', 'Aproximacion')
hold off;
legend;
title('Aproximacion s2 usando nodos arbitrarios');
xlabel('x');
ylabel('Valor de la funcion');

figure;
plot(xx, abs(f(xx) - yy))
disp(max(abs(f(xx) - yy)))
title('Error usando nodos arbitrarios');
xlabel('x');
ylabel('Error');

xn = linspace(0,1,m+1);
yn = f(xn);
sigmas = calculo_s2(xn, yn);
yy = arrayfun(@(x) s2(sigmas, x, xn, yn), xx);
figure;
plot(xx, f(xx),'', 'DisplayName', 'f(x)')
hold on;
plot(xx, yy, 'r--', 'DisplayName', 'Aproximacion')
hold off;
legend;
title('Aproximacion s2 usando nodos equisdistantes');
xlabel('x');
ylabel('Valor de la funcion');

figure;
plot(xx, abs(f(xx) - yy))
disp(max(abs(f(xx) - yy)))
title('Error usando nodos equisdistantes');
xlabel('x');
ylabel('Error');

function s2 = der_s2(sigmas, x, xn)
    %
    %   Función de la derivada explícita de s2
    %
    %   Inputs 
    %       sigmas: calculados con calculo_s2
    %       x: punto a evaluar
    %       xn: coordenadas x de los puntos
    %
    %   Outputs
    %       s2: valor de s2 en x 

    % Mismo principio de búsqueda
    nodes = find(xn>x);
    if x == 1 
        lower = length(xn)-1;
    else 
        lower = nodes(1)-1;
    end
    
    % Pero diferente fórmula, no ocupamos yn
    s1 = sigmas(lower);
    s2 = sigmas(lower+1);
    x1 = xn(lower);
    x2 = xn(lower + 1);
    s2 = s1 + ((s2 - s1)/((x2 - x1)))*(x-x1);
end

derf = @(x) -2.*pi.*sin(2.*pi.*x);
yy = arrayfun(@(x) der_s2(sigmas, x, xn), xx);
figure;
plot(xx, derf(xx),'', 'DisplayName', 'derf(x)')
hold on;
plot(xx, yy, 'r--', 'DisplayName', 'Aproximacion')
hold off;
legend;
title('Derivada s2 aproximada usando nodos equisdistantes');
xlabel('x');
ylabel('Valor de la funcion');

figure;
plot(xx, abs(derf(xx) - yy))
disp(max(abs(derf(xx) - yy)))
title('Error derivada s2 aprox usando nodos equisdistantes');
xlabel('x');
ylabel('Error');