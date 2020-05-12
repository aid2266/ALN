% Aidan Deaves

function [x,rho,res,iter] = jacobi(A,b,x0,nmax,prec)
% INPUT:        A, b matrices solve G.S
%               x0 vector inicial
%               nmax #max iteraciones, prec precisión
% OUTPUT:       Método de Jacobi: x vector sol, rho radio espectral               
%               res |Ax - b|_2, iter #iterations


% [0] Test matriz diagonal para valores pequeños (warning inversa!)
tol=1.0e-10;
D=diag(A);
if min(abs(D))<tol %check for zeros on the diagonal
    warning('error: values on the diagonal with abs. val < %e',tol)
end

% [0] Cálculos previos sobre los vectores INPUT
x = x0(:) % siempre devuleve vector columna! 
b = b(:); 
res0 = norm(b - A*x); 
if (res0 < prec) 
    warning('x0 such that |Ax0 - b| < prec!'); 
end

% [1] Aplicamos Jacobi
%     B = - inv(D) (L+U),     c = inv(D) * b
D = diag(D); 
L = tril(A, 1); % LT, diagonal de 0's
U = triu(A, 1); % UT, diagonal de 0's 
c = inv(D) * b;
B = - inv(D) * (L+U);   % equivalente: inv(L + D) * U

% Radio Espectral de B
rho = max(abs(eig(B))); 

% [2] Aplicamos método iterativo Gauss-Sidel
for iter = 1 : nmax
    x = B*x + c; 
    res = norm(b - A*x) / res0; 
    if (res < prec)
        return % hem acabat! 
    end
end

iter = -nmax; 
fprintf('No convergence in %d iterations', nmax); 
end
