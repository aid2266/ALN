
% aqui posarem la nostra resolucio del metode de GS modificat
function [Q,R] = gsm(A)
[m, n] = size(A); %% sabem nombre de columnes
Q = A; 
R = zeros(n); 

for k = 1:n
    R(k,k) = norm(Q(:,k), 2); 
    Q(:,k) = Q(:,k) / R(k,k) 
    kk = k+1; 
    R(k, kk:n) = transpose(Q(:, k)) * Q(:, kk:n); 
    Q(:, kk:n) = Q(:, kk:n) - Q(:, k)*R(k,kk:n); 
end

end