x = linspace(0,1,100);
y = linspace(0,2,200);

[xx,yy] = meshgrid(x,y);

f = @(x,y) sin(2.*pi.*(x+y)).*sin(pi.*(x-y));
surf(xx,yy,f(xx,yy))
view([-45 45])
shading interp
xlabel('xx')
ylabel('yy')
zlabel('f(x,y)')
title('Grafica de superficie de f(x, y)')

zz = f(xx,yy);
[U,S,V] = svd(zz);
val = diag(S);
aprox = cell(4,1);
aprox{1} = val(1)*U(:,1)*V(:,1)';
for i = 2:4
    aprox{i} = aprox{i-1} + val(i)*U(:,i)*V(:,i)';
end

figure(1)
for i = 1:4
    subplot(2,3,i);
    surf(xx, yy, aprox{i})
    shading interp
    xlabel('xx')
    ylabel('yy')
    zlabel('f(x,y)')
    title(sprintf('Aproximación de Rango %d', i))
end

subplot(2,3,5)
surf(xx,yy,f(xx,yy))
shading interp
xlabel('xx')
ylabel('yy')
zlabel('f(x,y)')
title('f(x, y)')

error = cell(4,1);
for i = 1:4
    error{i} = abs(f(xx, yy) - aprox{i});
    figure(2)
    subplot(2, 2, i)
    surf(xx, yy, error{i})
    shading flat
    xlabel('xx')
    ylabel('yy')
    zlabel('f(x,y)')
    title(sprintf('Error absoluto del Rango %d', i))
end

disp(val(1:4))
rank(zz)
   
disp(norm(cos(2.*pi.*x).*cos(pi.*x)) * norm(sin(2.*pi.*y).*sin(pi.*y)))
disp(norm(sin(2.*pi.*x).*sin(pi.*x)) * norm(cos(2.*pi.*y).*cos(pi.*y)))
disp(norm(sin(2.*pi.*x).*cos(pi.*x)) * norm(cos(2.*pi.*y).*sin(pi.*y)))
disp(norm(cos(2.*pi.*x).*sin(pi.*x)) * norm(sin(2.*pi.*y).*cos(pi.*y)))