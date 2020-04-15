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

void resol(MD& A, MD& A_Copy, VD& P, VD& b, VD& x, int N);
int find_max_pivot (const MD& A, int n, int k);
static bool abs_compare(double a, double b);
void write (MD& A, int n);
void writeV (VD& b, int n);

// ** FUNCION DESCOMPOSICION LU ** //

int lu(MD& A, VD& b, VD& x, int N, double tol) {
    /* IMPORTANTE NOTA AL LECTOR:
     este algoritmo para descomponer A en LU se realiza sobre la matriz A.
     La matriz descompuesta A' = L + U. Donde L es triangular inferior con 0's
     en la diagonal (añadimos 1's a posteriori) y U triangular superior.
     
     Se utiliza el pivotage parcial escalonado. Encontramos el pivote maximo
     mediante la funcion find_max_pivot que se encuentra al final del fichero.
    */
    
    MD L(N, VD(N));         // triangular inferior
    MD U(N, VD(N)) ;        // triangular superior
    VD P(N);                // vector de permutacion
    MD A_Copy(N, VD(N));    // copia de la matriz
    int perm_counter = 0;
    A_Copy = A;
    
    for (int i = 0; i < N; i++) P[i] = i; //init vect permutacion

    // ** FACTORIZACION LU sobre matriz A** //
    for (int i = 0; i < N; i++) {
        int imax = find_max_pivot(A, N, i); // pos pivote max, parcial escalonado, mirar funcion final del documento
        
        if (abs(A[imax][i]) < tol) return 0; //matriz es degenerada
        
        swap(A[i], A[imax]); // realizamos cambio de filas
        swap(P[i], P[imax]); // actualizamos vector perm.
        perm_counter++; // contador del numero de permutaciones

        for (int j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i]; // multiplicador
            for (int k = i + 1; k < N; k++){
                A[j][k] -= A[j][i] * A[i][k]; // algoritmo sobre A
            }
        }
    }
    
    // Hallamos L y U dentro de A descompuesta, recordamos A' = L + U
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j) {
                U[i][i] = A[i][i]; // pivotes en diagonal son de U
                L[i][i] = 1;       // imponemos diagonal de 1's
            }
            else if (i > j) L[i][j] = A[i][j];
            else U[i][j] = A[i][j];
        }
    }
    
    // ** CALCULO DEL DETERMINANTE ** //
    double det = 1.;
    for (int i = 0; i < N; i++) det *= U[i][i];
    
    // ** CALCULO DEL ERROR |PA - LU| norma 1 ** //
    MD LU(N, VD(N)); // inicializamos matriz producto L*U
    
    double norma_1 = 0.;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++){
                LU[i][j] += L[i][k] * U[k][j]; // calculamos producto LU
            }
            norma_1 = abs(A_Copy[P[i]][j] - LU[i][j]); // abs(PA - LU)
        }
    }
    
    // ** RESULTADOS ** //
    
    //cout << "Matriu A descomposada en LU" << endl;
    //write(LU, N);
    cout << "dim A: " << N << endl;
    cout << "error |PA - LU|_1: " << norma_1 << endl;
    cout << "Vector permutació: ";
    writeV(P, N);
    cout << "determinant: " << det << endl;
    
    resol(A, A_Copy, P, b, x, N); // llamamos funcion resolver
    
    cout << "funció lu retorna: ";
    if (perm_counter % 2 == 0) return 1; // existoso + #par de filas permutadas
    else return -1; // exitoso, #impar de filas permutadas
}


// ** CODIGO COMPLEMENTARIO A LA FUNCION LU ** //

// Escribe en pantalla una matriz A //
void write(MD& A, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl; cout << endl; cout << endl;
}

// Escribe en pantalla un vector b //
void writeV(VD& b, int n){
    cout << '[';
    for (int j = 0; j < n; j++){
            cout << b[j] << ' ';
        }
    cout << ']' << endl;
}


// Halla la posicion del pivote maximo, metodo parcial escalonado //
int find_max_pivot (const MD& A, int n, int k){
    
    /* IMPORTANTE NOTA AL LECTOR:
    La funcion max_element pertenece a la libreria algorithm que halla
    el elemento maximo de la fila de una matriz.
    Incluye como parametro la funcion abs_compare que compara los valores de
    cada fila en valor absoluto, extrayendo el maximo.
    */
    
    int pos_max = k; // indice de fila con pivote max
    double max_pivot = 0; // inicializamos a primer pivote
    
    for (int i = k; i < n; i++){
        double temp = *max_element(A[i].begin(), A[i].end(), abs_compare); // funcion que halla el elemento maximo de la fila,
        double elem = abs(A[i][k]/temp); // current element of column
        if (elem > max_pivot){
            max_pivot = elem;
            pos_max = i;
        }
    }
    return pos_max;
    
}

// ** Funcion que halla elemento maximo en valor absoluto ** //
static bool abs_compare(double a, double b){
    return (abs(double(a)) < abs(double(b)));
}

#endif


// AIDAN DEAVES - ALGEBRA LINEAL NUMERICA 2019 - 2020 FME
