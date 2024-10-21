% ejercicio 3
x = linspace(0,1,100);
y = linspace(0,2,200);

[xx,yy] = meshgrid(x,y);

f = @(x,y) sin(2.*pi.*(x+y)).*sin(pi.*(x-y));
surf(xx,yy,f(xx,yy))
view([-45 45])
shading interp

zz = f(xx,yy);
[U,S,V] = svd(zz);
val = diag(S);
aprox = cell(4,1);
aprox{1} = val(1)*U(:,1)*V(:,1)';
for i = 2:4
    aprox{i} = aprox{i-1} + val(i)*U(:,i)*V(:,i)';
end

for i = 1:4
    subplot(2,3,i);
    surf(xx, yy, aprox{i})
    shading interp
end

subplot(2,3,5)
surf(xx,yy,f(xx,yy))
shading interp

error = cell(4,1);
for i = 1:4
    error{i} = abs(f(xx, yy) - aprox{i});
    figure(2)
    subplot(2, 2, i)
    surf(xx, yy, error{i})
    shading flat
end

disp(val)
rank(zz)
    