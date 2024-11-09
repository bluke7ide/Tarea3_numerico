
%% Inciso a

m = 12;
x = ones(m, 1); % solución exacta
A = hilb(m);
b = A * x;

cond_A = cond(A);

%% Inciso b

[P, L, U] = lu(A);
y = L \ (P * b);
x_LU = U \ y;

error_LU = norm(x - x_LU, 2);

%% Inciso c 

[Q, R] = qr(A);
x_QR = R \ (Q' * b);

error_QR = norm(x - x_QR, 2);

%% Inciso d

L = chol(A, 'lower');
y = L \ b;
x_Chol = L' \ y;

error_Chol = norm(x - x_Chol, 2);

%% Inciso e

[U, S, V] = svd(A);
x_SVD = V * (S \ (U' * b));

error_SVD = norm(x - x_SVD, 2);

%% Inciso g

rango = 9;
A_nu = U(:, 1:rango) * S(1:rango, 1:rango) * V(:, 1:rango)';
x_aprox = A_nu \ b;

error_approx = norm(x - x_aprox, 2);
