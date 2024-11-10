
%% Inciso a

data = load("data/data.mat");

%% Inciso b

function fxc = funcion(x, c)
    % 
    %   Función que evalúa 
    %   f(x, c) = c1 * exp(-c2 * x) * sin(c3 * x + c4)
    %
    %   Inputs
    %       x: vector con las x
    %       c: vector con los coeficientes [c1, c2, c3, c4]
    %
    %   Output
    %       f: valores de f(x, c) evaluados en x
    
    fxc = c(1) * exp(-c(2) * x) .* sin(c(3) * x + c(4));
end

function Jf = jacobiano_f(x, c)
    %
    %   Función que calcula el Jacobiano de f respecto a los
    %   coeficientes c = [c1, c2, c3, c4].
    %
    %   Inputs
    %       x: vector con las x
    %       c: vector con los coeficientes [c1, c2, c3, c4]
    %
    %   Output
    %       J: matriz Jacobiana evaluada en x y c
    
    Jf = zeros(length(x), 4);

    Jf(:, 1) = exp(-c(2) * x) .* sin(c(3) * x + c(4));        
    Jf(:, 2) = -c(1) * x .* exp(-c(2) * x) .* sin(c(3) * x + c(4));    
    Jf(:, 3) = c(1) * exp(-c(2) * x) .* x .* cos(c(3) * x + c(4));
    Jf(:, 4) = c(1) * exp(-c(2) * x) .* cos(c(3) * x + c(4));    
end

function [c, iter, norma_delta_c, num_cond] = minimos_cuadrados(x, y, c0)
    %
    %   Función que ejecuta el algoritmo para minimizar la diferencia
    %   entre la función f(x, c) y los valores observados y.
    %
    %   Inputs
    %       x: vector con las x
    %       y: vector con las y
    %       c0: vector con los coeficientes [c1, c2, c3, c4]
    %
    %   Outputs
    %       c: vector de coeficientes ajustados
    %       iter: número de iteraciones requeridas
    %       norma_delta_c: valor de la norma de delta c por iteración
    %       num_cond: número de condición por itereción
    
    
    tol = 1e-8;
    c = c0;
    iter = 0;
    delta_c = inf;
    num_cond = [];
    
    while norm(delta_c, inf) > tol
       
        dif = y - funcion(x, c);
        J = jacobiano_f(x, c);
        
        delta_c = (J' * J) \ (J' * dif);
        c = c + delta_c;

        iter = iter + 1; 
        norma_delta_c(iter) = norm(delta_c, inf);
        num_cond(iter) = cond(J);
    end
end

%% Inciso c

x = data.x;
y = data.y;
c0 = [1.1; 0.4; 2.1; 0.2];

% Ejecución del algoritmo
[c, iter, norma_delta_c, num_cond] = minimos_cuadrados(x, y, c0);

disp("Número de iteraciones:");
disp(iter);
disp("Coeficientes ajustados:");
disp(c);

% Función exacta
c_exact = [1; 0.5; 2; 0];
xx = linspace(0, 10);
yy = funcion(xx, c_exact);

% Gráficas
figure;
hold on;
plot(x, y, 'b.', 'DisplayName', 'Datos (x_i, y_i)', 'MarkerSize', 10);
plot(xx, yy, 'k-', 'DisplayName', 'Función Exacta');
plot(sort(x), funcion(sort(x), c), 'r--', 'DisplayName', 'Aproximación', 'LineWidth', 1.5);
legend;
xlabel('x');
ylabel('y');
title('Problema de mínimos cuadrados');
hold off;

figure;
plot(1:iter, norma_delta_c, '.-', 'MarkerSize', 15);
xlabel('Iteración k');
ylabel('||\Deltac_k||_\infty');
title('||\Deltac_k||_\infty en función de k');
grid on;


%% Inciso d

figure;
plot(1:iter, num_cond, '.-', 'MarkerSize', 15);
xlabel('Iteración k');
ylabel('\kappa(J(c_k))');
title('Número de condición en función de k');
grid on;

%% Inciso e

c0 = [1; 1; 1; 1];
[c, iter] = minimos_cuadrados(x, y, c0);

disp("Número de iteraciones:");
disp(iter);
disp("Coeficientes ajustados:");
disp(c);


