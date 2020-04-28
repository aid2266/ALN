format long; % això fa que la sortida sigui amb més decimals
clear vars all;
% llegim d'un fitxer
A = xlsread('world_population_1800.xls');
year=A(:,1); % year a la primera columna
popul=A(:,2); % population segona columna
n=length(year);  % longitud del vector year

% fem un primer dibuix
figure(1)
scatter(year,popul,'r','filled'); % dibuixa punts
hold on; % això permet continuar dibuixant sobre aquesta figura
         % més endavant. Per acabar: hold off

% volem resoldre p = A e^By
% prenem logs --> z = log(p), z = log(A) + By
% per tant volem resoldre z = alpha + By
z = log(popul); %log nep. 

%Fem descomposició QR 
%%R = zeros(2,2);
%Q = zeros(12, 2); 

%Construim M
M = [ones(n,1) , year]; 
Mcopia = M; 

%primer pas GS Mod
R(1,1) = norm(M(:,1),2); 
M(:, 1) = M(:,1)/R(1,1); 

% pas 2
R(1,2) = dot(M(:,1),M(:,2)); 
M(:, 2) = M(:,2) - R(1,2)*M(:,1); %segona columna normalitzada


%normalitzar columna 2 
R(2,2) = norm(M(:,2),2);
M(:,2) = M(:,2)/R(2,2); 
Q = M; 
Q
R
%comprovem que q es ortogonal
transpose(Q)*Q
% qTq = Id

% Resolem les equacions normals
sol = R\(transpose(Q)*z) % vector solucio de la forma [alpha + B]
alpha = sol(1)
B = sol(2)

%dibuixar sobre figure(1) la grafica A*exp(B * year)
% comparem punts i solucio

figure(2) 
scatter(year, z, 'r', 'filled'); 
hold on; 
yy = ((year(1) - 0.01) : 0.01 : (year(n) + 0.01)); %vector 
zaprox = alpha + B*yy; 
plot(yy, zaprox, '-')




