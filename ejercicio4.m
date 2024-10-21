% ejercicio 4
function c = calculo_aprox(f, n)
    M = zeros(n);
    for i = 1:n
        for j = 1:n
            M(i,j) = 1/(i+j+1);
        end
    end

    
end

f = @(x) cos(4.*pi.*x);
calculo_aprox(f, 10)
