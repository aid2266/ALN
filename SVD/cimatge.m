clearvars
close all
clc

A = imread('koshka.jpg'); % aquesta es la imatge en RGB, size 500 x 500 x 3
A = rgb2gray(A); 
A = double(A); 
r = rank(A) %tinc r valors singulars diferents de la tolerancia

k = [5, 10, 25, 50, 75, 100]; %nombre de valors singulars
[m n] = size(A); 
sizeA = m * n; %nombre de components a la imatge

% fem la descomposicio SVD
[U, S, V] = svd(A); 

for i = k % nombre de components de k
    X = U(:, 1:i) * S(1:i, 1:i) * V(:, 1:i)'; 
    figure() %obro finesta grafica
    imshow(uint8(X)) %passem X de double a enter, 
    display('press intro')
    %title ('#SV: %d', num2str(i), ' of ', num2str(r)); 
    pause
end

imshow(uint8(A))
pause; 
