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
    x(end, 1) = 1;
    x(end, end) = -1;
    y(end) = 0;

    sigmas = x\y;
    
end

m = 15;
xn = sort(rand(m+1,1));
xn(1) = 0;
xn(end) = 1;
yn = f(xn);
f = @(x) cos(2.*pi.*x);
sigmas = calculo_s2(xn, yn);

function s2 = s2(sigmas, x, xn, yn)
    nodes = find(xn>x);
    if x == 1 
        lower = length(xn)-1;
    else 
        lower = nodes(1)-1;
    end
    
    s1 = sigmas(lower);
    s2 = sigmas(lower+1);
    y1 = yn(lower);
    x1 = xn(lower);
    x2 = xn(lower + 1);
    s2 = y1 + s1.*(x - x1) + ((s2 - s1)/(2 * (x2 - x1)))*(x-x1)^2;
end

% s2(sigmas, 1, xn, yn)
xx = linspace(0,1);
yy = arrayfun(@(x) s2(sigmas, x, xn, yn), xx);
plot(xx, f(xx))
hold on;
plot(xx, yy)
hold off;

figure(2)
plot(xx, abs(f(xx) - yy))
disp(max(abs(f(xx) - yy)))

xn = linspace(0,1,m+1);
yn = f(xn);
sigmas = calculo_s2(xn, yn);
yy = arrayfun(@(x) s2(sigmas, x, xn, yn), xx);
figure(3)
plot(xx, f(xx))
hold on;
plot(xx, yy)
hold off;

figure(4)
plot(xx, abs(f(xx) - yy))
disp(max(abs(f(xx) - yy)))

function s2 = der_s2(sigmas, x, xn)
    nodes = find(xn>x);
    if x == 1 
        lower = length(xn)-1;
    else 
        lower = nodes(1)-1;
    end
    
    s1 = sigmas(lower);
    s2 = sigmas(lower+1);
    x1 = xn(lower);
    x2 = xn(lower + 1);
    s2 = s1 + ((s2 - s1)/((x2 - x1)))*(x-x1);
end

derf = @(x) -2.*pi.*sin(2.*pi.*x);
yy = arrayfun(@(x) der_s2(sigmas, x, xn), xx);
figure(5)
plot(xx, derf(xx))
hold on;
plot(xx, yy)
hold off;

figure(6)
plot(xx, abs(derf(xx) - yy))
disp(max(abs(derf(xx) - yy)))