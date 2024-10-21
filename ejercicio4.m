% ejercicio 4
function c = calculo_aprox(f, n, inf, sup)
    M = zeros(n);
    b = zeros(n);
    for i = 1:n
        b(i) = integral(@(x) f*x^(i-1), inf, sup);
        for j = 1:n
            M(i,j) = 1/(i+j+1);
        end
    end
    
    % M tiene rango completo, procedemos por pseudo-inversa
    inversa = (M'\M)*M';
    c = inversa*b;

end

f = @(x) cos(4.*pi.*x);
calculo_aprox(f, 10)
