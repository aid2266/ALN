#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
//#include "lu.h"
//#include "resol.h"

using namespace std;

const double tol = 0.0000000000001;

typedef vector<double> VD;
typedef vector <VD> MD;
int LUPDecompose(MD& A, VD& b, VD& x, int n, double tol);
void LUPSolve(MD& A, MD& A_Copy, VD& P, VD& b, VD& x, int N);
int find_max_pivot (const MD& A, int n, int k);
static bool abs_compare(double a, double b);
void write (MD& A, int n);
void writeV (VD& b, int n);


int main(){
    
    ifstream inFile;
    inFile.open("/Users/aidandeaves/Documents/ALN/ALN_P1/ALN_P2/MAT/M02.DAT"); // change this file!!
    
    // check for any errors
    if (inFile.fail()){
        cerr << "Error Llegint Fitxer" << endl;
        return 0;
    }
    
    int n; // dim del sistema cuadrado
    int v; // numero de comoponentes no nulos de la matriz A
    int k; // ''     ''     ''       ''   ''  de la matriz b
    
    inFile >> n >> v;
    MD A(n, VD(n));     // creamos matriz A
    VD b(n);            // vector b
    VD x(n);            // vector solucion
    
    for (int counter = 0; counter < v; counter++){
        int i, j;
        double elem;
        inFile >> i >> j >> elem;
        A[i][j] = elem;
    }
    
    inFile >> k;
    
    for (int counter = 0; counter < k; counter++){
        int i;
        double elem;
        inFile >> i >> elem;
        b[i] = elem;
    }
    
    cout << LUPDecompose(A, b, x, n, tol) << endl; // llamamos funcion descomponer
    
    // ** output.txt ** //
    // tenemos la solucion pasada por referencia
    
}


int LUPDecompose(MD& A, VD& b, VD& x, int N, double Tol) {

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

    // ** FACTORIZACION LU guardado en la matriz A** //
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
    cout << "Matriz A descompuesta en LU" << endl;
    write(LU, N);
    cout << "dim A: " << N << endl;
    cout << "error |PA - LU|_1: " << norma_1 << endl;
    cout << "Vector permutacion: ";
    writeV(P, N);
    cout << "determinante: " << det << endl;
    
    LUPSolve(A, A_Copy, P, b, x, N); // llamamos funcion resolver
    
    if (perm_counter % 2 == 0) return 1; // existoso + #par de filas permutadas
    else return -1; // exitoso, #impar de filas permutadas
}


// ** RESOLVEMOS LUx = b ** //

void LUPSolve(MD& A, MD& A_Copy, VD& P, VD& b, VD& x, int N) {
    
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
        x[i] = x[i] / A[i][i];
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
    
    
    cout << "solucio del sistema" << endl;
    writeV(x, N);
    
}


// ** FUNCIONES COMPLEMENTARIAS ** //

// Escribe la matriz //
void write(MD& A, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl; cout << endl; cout << endl;
}

// Escribe un vector //
void writeV(VD& b, int n){
    cout << '[';
    for (int j = 0; j < n; j++){
            cout << b[j] << ' ';
        }
    cout << ']' << endl;
}


// Encuentra pivote maximo //
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



