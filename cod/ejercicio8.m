% ejercicio 8

function sigmas = calculo_s2(xn, yn)
    m = length(xn);
    x = zeros(m);
    y = zeros(m,1);
    for i = 1:(m-1)
        x(i,i) = (xn(i+1)-xn(i))/2;
        x(i,i+1) = x(i,i);
        y(i) = yn(i+1) - yn(i);
    end
    x(end, end) = 1;
    y(end) = 0;

    sigmas = x\y;
    
end

m = 15;
xn = sort(rand(m+1,1));
xn(1) = 0;
xn(end) = 1;
f = @(x) cos(2.*pi.*x);
sigmas = calculo_s2(xn, f(xn));