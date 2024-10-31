% ejercicio 4
function c = calculo_aprox(f, n, inf, sup)
    M = zeros(n);
    b = zeros(n,1);
    for i = 1:n
        b(i) = integral(@(x) f(x).*x.^(i-1), inf, sup);
        for j = 1:n
            M(i,j) = 1/(i+j+1);
        end
    end
    
    % M tiene rango completo, procedemos por pseudo-inversa
    inversa = (M'\M)*M';
    c = inversa*b;

end

f = @(x) cos(4.*pi.*x);
c = calculo_aprox(f, 10, 0, 1);
p = @(x) polyval(c(end:-1:1), x);
integral(@(x) f(x) - p(x), 0, 1)

xx = linspace(0,1);
plot(xx, f(xx), '-')
hold on
plot(xx, p(xx), 'b-')
hold off



