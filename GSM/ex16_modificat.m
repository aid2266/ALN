format long; % number of digits set to long

%% construimos la matriz con la que vamos a operar
x = (0.25 : 0.25 : 1.75);
y = [0.4 ;0.5; 0.9; 1.28; 1.6; 1.66; 2.02];

%hacemos un dibujo de los datos
figure(1)
scatter(x, y, 'r', 'filled'); 
hold on; 

m = length(x); 
deg = 6; %% m'especifica el grau del polinomi
n = deg + 1; 

A = vander(x); 
A = A(:, m-n+1:m);

%% buscamos una solucion del tipo P(x) = a + a_1x + a_2x^2

%Inicializamos y hacemos copias de la matriz
A_Copy = A; 
Q = A;
R = zeros(n); 

% Llamamos a la funcion gsm que nos calcula Gram-Schmidt Modificado
[Q,R] = gsm(A); 

error = norm(Q'*Q - eye(n)) %test de la identitat

%% resolem el sistema
sol = R\(transpose(Q)*y); 

%scatter plot 2
h = 0.001; 
xx = [0:h:2];
yy = polyval(sol, xx); 
plot(xx, yy, '-r'); %% esta esta mal

%residu del sistema
res1 = norm(A_Copy*sol - y, 2)
