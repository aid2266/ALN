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
int lu(MD& A, VD& b, int n, double tol);
int LUPDecompose(MD& A, VD& b, int n, double tol);
void LUPSolve(MD& A, MD& A_Copy, VD& P, VD& b, int N);
int copyLUP(MD& A, VD& b, int n, double tol);
void resol(MD& L, MD& U, VD& b, VD& perm, MD& A, int n);
MD multiplyMatrix(MD& A, MD& B, int n);
int find_max_pivot (const MD& A, int n, int k);
static bool abs_compare(double a, double b);
void write (MD& A, int n);
void writeV (VD& b, int n);

int main(){
    
    ifstream inFile;
    inFile.open("/Users/aidandeaves/Documents/ALN/ALN_P1/ALN_P2/MAT/M01.DAT"); // change this file!!
    
    // check for any errors
    if (inFile.fail()){
        cerr << "Error Obrint Fitxer" << endl;
        return 0;
    }
    
    int n; // dim del sistema cuadrado
    int v; // numero de comoponentes no nulos de la matriz A
    int k; // ''     ''     ''       ''   ''  de la matriz b
    
    inFile >> n >> v;
    MD A(n, VD(n));     // creamos las matrices correspondientes
    VD b(n);
    
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
    
    cout << LUPDecompose(A, b, n, tol) << endl; // call the function
    
    // ** una vez realizada la descomposicion LU ** //
    cout << "dimension del sistema: " << n << endl;
    cout << "estimacion del error de descomposicion PA = LU: " << 'e' << endl;
    cout << "vector de permutaciones de P: " << 'p' << endl;
    cout << "determinant de la matriu A: " << 'd' << endl;
    cout << "estimacio del error del sistema Ax = b, amb norma infinita: " << 'i' << endl;
    
}


int LUPDecompose(MD& A, VD& b,  int N, double Tol) {

    /* DISCLAIMER:
     este algoritmo para descomponer A en LU se realiza sobre la matriz A.
     La matriz descompuesta A' = L + U. Donde L es una mat-tri-inf con 0's
     en la diagonal (añadimos 1's a posteriori) y U una mat-tri-sup.
    */
    
    int imax;
    int perm_counter = 0;
    MD L(N, VD(N));     // triangular inferior
    MD U(N, VD(N)) ;    // triangular superior
    VD P(N);               // vector de permutacion
    MD A_Copy(N, VD(N));   // copia de la matriz
    A_Copy = A;
    
    for (int i = 0; i < N; i++) P[i] = i; //init vect permutacion unitario

    // ** FACTORIZACION LU guardado en la matriz A** //
    for (int i = 0; i < N; i++) {
        int imax = find_max_pivot(A, N, i); // pos pivote max, parcial escalonado
        if (abs(A[imax][i]) < tol) return 0; //matriz es degenerada
        swap(A[i], A[imax]);
        swap(P[i], P[imax]);
        perm_counter++; // contador del numero de permutaciones

        for (int j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i]; // calculo del multiplicador

            for (int k = i + 1; k < N; k++){
                A[j][k] -= A[j][i] * A[i][k];
            }
        }
    }
    
    // Hallamos L y U dentro de A descompuesta, recordamos A' = L + U
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j) {
                U[i][i] = A[i][i];
                L[i][i] = 1; // metodo Doolitle
            }
            else if (i > j) L[i][j] = A[i][j];
            else U[i][j] = A[i][j];
        }
    }
    
    // ** CALCULO DEL DETERMINANTE ** //
    double det = 1.;
    for (int i = 0; i < N; i++) det *= U[i][i];
    
    // ** CALCULO DEL ERROR |PA - LU| norma 1 ** //
    MD LU(N, VD(N)); // matriz producto L*U
    
    double norma_1 = 0.;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++){
                LU[i][j] += L[i][k] * U[k][j]; // calculamos LU
            }
            norma_1 = abs(A_Copy[P[i]][j] - LU[i][j]); // PA - LU
        }
    }
    cout << "Matriz A descompuesta en LU" << endl;
    write(LU, N);
    cout << "Vector permutacion: " << '[';
    writeV(P, N);
    cout << ']' << endl;
    cout << "dim A: " << N << endl;
    cout << "determinante: " << det << endl;
    cout << "error |PA - LU|_1: " << norma_1 << endl;
    
    LUPSolve(A, A_Copy, P, b, N);
    
    if (perm_counter % 2 == 0) return 1; // numero par de ops
    else return -1;              // numero impar de ops
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
void LUPSolve(MD& A, MD& A_Copy, VD& P, VD& b, int N) {
    
    VD y(N);
    VD x(N); // creamos el vector solucion
    // tenemos que resolver LUx = b
    
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
    double norma_inf = 0.;
    
    for (int i = 0; i < N; i++){
        double Ax = 0.;
        for (int j = 0; j < N; j++){
            Ax += A_Copy[i][j]*x[j]; // valor de b'
        }
        norma_1 += Ax - b[i]; //
    }
    
    cout << "norma 1 error |Ax* - b|: " << norma_1 << endl;
    
    
    cout << "solucion del sistema" << endl;
    for (int i = 0; i < N; i++) cout << x[i] << ' ';
    cout << endl;
    
}

/* INPUT: A,P filled in LUPDecompose; N - dimension
 * OUTPUT: IA is the inverse of the initial matrix
 */
void LUPInvert(double **A, int *P, int N, double **IA) {
  
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            if (P[i] == j)
                IA[i][j] = 1.0;
            else
                IA[i][j] = 0.0;

            for (int k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }

        for (int i = N - 1; i >= 0; i--) {
            for (int k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];

            IA[i][j] = IA[i][j] / A[i][i];
        }
    }
}

/* INPUT: A,P filled in LUPDecompose; N - dimension.
 * OUTPUT: Function returns the determinant of the initial matrix
 */
double LUPDeterminant(double **A, int *P, int N) {

    double det = A[0][0];

    for (int i = 1; i < N; i++)
        det *= A[i][i];

    if ((P[N] - N) % 2 == 0)
        return det;
    else
        return -det;
    
    cout << "determinante del sistema " << det << endl;
}

// escribe matriz
void write(MD& A, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl; cout << endl; cout << endl;
}

// escribe vector
void writeV(VD& b, int n){
    cout << '[';
    for (int j = 0; j < n; j++){
            cout << b[j] << ' ';
        }
    cout << ']' << endl;
    
}

int find_max_pivot (const MD& A, int n, int k){
    int pos_max = k; // index of the row with max pivot
    double max_pivot = 0; // assume it is initial pivot
    
    for (int i = k; i < n; i++){
        double temp = *max_element(A[i].begin(), A[i].end(), abs_compare); // max element of row
        double elem = abs(A[i][k]/temp); // current element of column
        if (elem > max_pivot){
            max_pivot = elem;
            pos_max = i;
        }
    }
    return pos_max;
    
}

// ** Function that compares absolute values results ** //
static bool abs_compare(double a, double b){
    return (abs(double(a)) < abs(double(b)));
}



