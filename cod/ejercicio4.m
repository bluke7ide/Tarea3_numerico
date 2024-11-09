% ejercicio 4
function c = calculo_aprox(f, n, inf, sup)
    M = zeros(n);
    b = zeros(n,1);
    for i = 1:n
        b(i) = integral(@(x) f(x).*x.^(i-1), inf, sup);
        for j = 1:n
            M(i,j) = 1/(i+j-1);
        end
    end
    
    c = M\b;

end

f = @(x) cos(4.*pi.*x);
c = calculo_aprox(f, 10, 0, 1);
p = @(x) polyval(c(end:-1:1), x);
integral(@(x) (f(x) - p(x)).^2, 0, 1)

xx = linspace(0,1); 
plot(xx, f(xx), '-', 'DisplayName', 'f(x) = cos(4\pi x)')
hold on
plot(xx, p(xx), 'r--', 'DisplayName', 'Aproximación p(x)')
hold off
legend;
title('Aproximación de f(x) = cos(4\pi x) usando derivadas');
xlabel('x');
ylabel('Valor de la función');

function p = aprox_leg(n, f)
    pn = cell(n, 1);
    pn{1} = @(x) 1; 
    pn{2} = @(x) x; 
    
    for k = 2:(n-1)
        pn{k+1} = @(x) ((2*k-1)*x.*pn{k}(x) - (k-1)*pn{k-1}(x)) / (k);
    end
  
    betas = zeros(n, 1);
    for k = 1:n
        betas(k) = integral(@(x) f(x) .* pn{k}(2*x-1), 0, 1) / ...
                   integral(@(x) pn{k}(2*x-1).^2, 0, 1, 'ArrayValued', true);
    end
    p = @(x) sum(cell2mat(arrayfun(@(k) betas(k) .* pn{k}(2*x-1), ...
        1:n, 'UniformOutput', false)));
end

n = 10;
p = aprox_leg(n, f);

figure(2)
xx = linspace(0, 1, 100); 
plot(xx, f(xx), '-', 'DisplayName', 'f(x) = cos(4\pi x)'); 
hold on;
plot(xx, arrayfun(p, xx), 'r--', 'DisplayName', 'Aproximación p(x)');
hold off;
legend;
title('Aproximación de f(x) = cos(4\pi x) usando polinomios de Legendre');
xlabel('x');
ylabel('Valor de la función');

norma_error = integral(@(x) (f(x) - arrayfun(p, x)).^2, 0, 1);
disp(['Norma del error: ', num2str(norma_error)]);