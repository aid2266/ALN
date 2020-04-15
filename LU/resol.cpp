// AIDAN DEAVES - ALGEBRA LINEAL NUMERICA 2019 - 2020 FME


#ifndef MAIN
#define MAIN

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;

typedef vector<double> VD;
typedef vector <VD> MD;

void writeV(VD& b, int n);

void resol(MD& A, MD& A_Copy, VD& P, VD& b, VD& x, int N) {
    
    VD y(N); // Ly = b
    
    // ** RESOLVER Ly = b, forward-substitution ** //
    for (int i = 0; i < N; i++) {
        y[i] = b[P[i]]; // solucion permutada
        for (int k = 0; k < i; k++)
            y[i] -= A[i][k] * y[k];
    }
    
    // ** RESOLVER Ux = y, back-substitution ** //
    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int k = i + 1; k < N; k++) x[i] -= A[i][k] * x[k];
        x[i] = x[i] / A[i][i]; // pasamos por referencia a main este resultado
    }

    // ** CALCULO DE ERRORES |Ax* - b| ** //
    double norma_1 = 0.;
    double norma_2 = 0.;
    double max_norma = 0;
    
    for (int i = 0; i < N; i++){
        double Ax = 0.; // calculo de Ax (matriz original por vector solucion)
        for (int j = 0; j < N; j++){
            Ax += A_Copy[P[i]][j]*x[j]; // utilizamos matriz original permutada
        }
        norma_1 += abs(Ax - b[P[i]]);
        norma_2 += norma_1 * norma_1;
        
        double temp = abs(Ax - b[P[i]]);
        if (temp > max_norma) max_norma = temp; // hallamos fila maxima
    }
    
    cout << "norma 1 error |Ax* - b|: " << norma_1 << endl;
    cout << "norma 2 error |Ax* - b|: " << sqrt(norma_2) << endl;
    cout << "norma infinit error |Ax* - b|: " << max_norma << endl;
}

#endif

// AIDAN DEAVES - ALGEBRA LINEAL NUMERICA 2019 - 2020 FME
