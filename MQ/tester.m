
clearvars all
% file with info for polminquad

x = [1800 1850 1875 1900 1920 1930 1940 1950 1960 1970 1980 1990]
y = [900 1200 1325 1625 1813 1987 2213 2516 3019 3693 4450 5333]
grau = 2
plt = 1

[coefs, norm2Res] = polminquad(x, y, grau, plt)
