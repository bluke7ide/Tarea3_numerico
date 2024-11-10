
%% Inciso a

function [Q, R] = qr1(A)
    %
    %   Función que calcula la factorización QR de una matriz A
    %   utilizando el método de Gram-Schmidt inestable.
    %
    %   Inputs
    %       A: Matriz de tamaño mxn
    %
    %   Outputs
    %       Q: Matriz unitaria de tamaño mxn
    %       R: Matriz triangular superior de tamaño nxn
    %

    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);

    for j = 1:n
        vj = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            vj = vj - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(vj);
        Q(:, j) = vj / R(j, j);
    end
end


%% Inciso b

function [Q, R] = qr2(A)
    %
    %   Función que calcula la factorización QR de una matriz A
    %   utilizando el método de Gram-Schmidt estable.
    %
    %   Inputs
    %       A: Matriz de tamaño mxn
    %
    %   Outputs
    %       Q: Matriz unitaria de tamaño mxn
    %       R: Matriz triangular superior de tamaño nxn
    %

    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    V = A;

    for i = 1:n
        R(i, i) = norm(V(:, i));
        Q(:, i) = V(:, i) / R(i, i);
        for j = i+1:n
            R(i, j) = Q(:, i)' * V(:, j);
            V(:, j) = V(:, j) - R(i, j) * Q(:, i);
        end
    end
end


%% Inciso c

function [Q, R] = qr3(A)
    %
    %   Función que calcula la factorización QR de una matriz A
    %   utilizando la triangularización de Householder.
    %
    %   Inputs
    %       A: Matriz de tamaño mxn
    %
    %   Outputs
    %       Q: Matriz unitaria de tamaño mxm
    %       R: Matriz triangular superior de tamaño mxn
    %

    [m, n] = size(A);
    R = A;
    Q = eye(m);

    for j = 1:n
        x = R(j:m, j);
        e1 = zeros(length(x), 1);
        e1(1) = 1;
        vj = sign(x(1)) * norm(x) * e1 + x;
        vj = vj / norm(vj);
        R(j:m, j:n) = R(j:m, j:n) - 2 * vj * (vj' * R(j:m, j:n));
        Q(j:m, :) = Q(j:m, :) - 2 * vj * (vj' * Q(j:m, :));
    end

    Q = Q';
end


%% Inciso d

m = 20; 
A_random = rand(m);

[Q1, R1] = qr1(A_random);
[Q2, R2] = qr2(A_random);
[Q3, R3] = qr3(A_random);

error_random = zeros(3,2);

error_random(1,1) = norm(A_random - Q1 * R1, 2);
error_random(1,2) = norm(Q1 * Q1' - eye(m), 2);

error_random(2,1) = norm(A_random - Q2 * R2, 2);
error_random(2,2) = norm(Q2 * Q2' - eye(m), 2);

error_random(3,1) = norm(A_random - Q3 * R3, 2);
error_random(3,2) = norm(Q3 * Q3' - eye(m), 2);


%% Inciso e

A_hilbert = hilb(m);

[Q1, R1] = qr1(A_hilbert);
[Q2, R2] = qr2(A_hilbert);
[Q3, R3] = qr3(A_hilbert);

error_hilbert = zeros(4,2);

error_hilbert(1,1)  = norm(A_hilbert - Q1 * R1, 2);
error_hilbert(1,2) = norm(Q1 * Q1' - eye(m), 2);

error_hilbert(1,2) = norm(A_hilbert - Q2 * R2, 2);
error_hilbert(2,2)  = norm(Q2 * Q2' - eye(m), 2);

error_hilbert(3,1) = norm(A_hilbert - Q3 * R3, 2);
error_hilbert(3,2) = norm(Q3 * Q3' - eye(m), 2);

% Comparación con la factorización QR de MATLAB
[Q4, R4] = qr(A_hilbert);

error_hilbert(4,1) = norm(A_hilbert - Q4 * R4, 2);
error_hilbert(4,2) = norm(Q4 * Q4' - eye(m), 2);

