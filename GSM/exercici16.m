format long; % number of digits set to long

%construimos la matriz con la que vamos a operar
v = (0.25 : 0.25 : 1.75);
b = [0.4 ;0.5; 0.9; 1.28; 1.6; 1.66; 2.02];

%hacemos un dibujo de los datos
figure(1)
scatter(v, b, 'r', 'filled'); 
hold on; 

%creamos la matriz A mediante potencias de las columnas
A = fliplr(vander(v)); 
A = A(:, 1:3);
n = 3; 

% buscamos una solucion del tipo P(x) = a + a_1x + a_2x^2

%Hacemos copia de la matriz
A_Copy = A; 

% haremos QR sobre la matriz A
% Paso 1
R(1,1) = norm(A(:,1),2); 
A(:,1) = A(:,1) / R(1,1); %normalizamos la columna 1

% podem posar
R(1,2:n) = A(:,1)' * A(:,2:n);  %dot no funciona
A(:,2:n) = A(:,2:n) - A(:,1)*R(1, 2:n); 

% Paso 2
R(2,2) = norm(A(:,2),2); 
A(:,2) = A(:,2) / R(2,2); 
R(2,3) = dot(A(:,2),A(:,3)); 
A(:,3) = A(:,3) - A(:,2)*R(2,3); %ortog. columna 3 respecto 2

% Paso 3
R(3,3) = norm(A(:,3),2); 
A(:,3) = A(:,3) / R(3,3); 

Q = A;
R; 
transpose(Q)*Q %test de la identitat

%resolem sistema 
sol = R\(transpose(Q)*b); 
sol_flip = [sol(3) sol(2) sol(1)] % le damos la vuelta

%scatter plot 2
h = 0.001; 
x = [0:h:2];
yy = polyval(sol_flip, x); 
plot(x, yy, '-r'); %% esta esta mal

%residu del sistema
res1 = norm(A_Copy*sol - b, 2)

%% evaluar polinomio en un intervalo con polyval
p = polyval(fliplr(sol), v); 
res2 = norm(p - b, 2) % no em funciona! 


