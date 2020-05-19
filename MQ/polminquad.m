% Aidan Deaves

function [coefs, norm2Res] = polminquad(x, y, grau, plt)
%  INPUT:      x - vector with x-coordinate data of Table 1
%              y - vector with y-coordinate data of Table 1
%              grau - degree of fitting polynomial by least squares
%              plt - graph of Table 1 data + smooth plot of fitting polynomial

% OUTPUT:      coefs - vector with least squares solution coefficients
%              norm2Res - residue of |Aa - y| using the 2-norm

%% [1] Creating matrix A of size m x n to QR decompose 

x = x(:); % change x and y to column vectors
y = y(:); 
m = length(x); 
n = grau + 1; 
A = vander(x); % Vandermonde matrix
A = A(:, m-n+1:m); % imposes m x n Vandermonde Matrix size
A_Copy = A; %copy of matrix for later use

%% [2] QR decompose matrix A using Modified Gram-Schmidt process
[Q, R] = GSM(A); % returns Q orthog. m x n matrix 
                 % R upper triangular n x n matrix
                 
%% [3] Test for Q's orthogonality |Q'Q - I| inf-norm
Q_res = norm(Q'* Q - eye(n), inf)
                 
%% [4] Solve linear system Ra = Q'y, checks if R is singular
[a] = solveQR(Q, R, y); 
coefs = a; % returns value of solution coefficients 
norm2Res = norm(A_Copy*a - y, 2); % least squares error

%% [5] Evaluate fitting polynomial + plot
h = 0.01; %interval mimicks 'smoothness'
xx = [x(1)-10*h:h:x(end)+10*h]; %space for graphing
yy = polyval(a, xx); %evaluate fitting with least squares solution

if plt == 1 % plot graph if true
    figure(1); 
    title('Table 1 data') 
    xlabel('x')
    ylabel('y')
    scatter(x, y, 'r', 'filled') %plots x and y  
    hold on; 
    plot(xx, yy); % plots fitting polynomial
end

end

%% Modified Gram-Schmidt Function
function [Q R] = GSM(A)
%  INPUT:     Matrix A of size m x n to QR decompose
%  OUTPUT:    Q, R matrices such that A = QR
[m n] = size(A);
Q = A;
R = zeros(n); % init matrix of zeros size n x n 

for k = 1:n
    R(k,k) = norm(Q(:,k), 2); %norm of kth column
    Q(:,k) = Q(:,k) / R(k,k); %normalises kth column
    j = k + 1; 
    R(k,j:end) = Q(:,k)' * Q(:,j:end);
    Q(:,j:end) = Q(:,j:end) - Q(:,k)*R(k,j:end); %orthogonalise k+1 columns
end  
end


%% Solves upper triangular system Ra = Q'y, checks if R is singular
function [a] = solveQR(Q, R, y)
% INPUT:    Matrices Q, R from mod. Gram-Schmidt process
%           y-axis vector
% OUTPUT:   least squares solution vector 'a'
%           EXTRA: checks if R is upper triangular

% Checks if R is singular, ie det is smaller than TOL
TOL = 1.e-12; 
detR = prod(diag(R)); 
if (detR < TOL) 
    error('R is singular, det R = %e (TOL %e)',det, TOL)
end

% Solve system
a = R\(transpose(Q) * y); 

end 
