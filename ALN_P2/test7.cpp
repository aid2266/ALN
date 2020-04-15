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
void LUPSolve(MD& A, VD& P, VD& b, int N);
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

    int i, j, k, imax;
    MD L(N, VD(N));     // triangular inferior
    MD U(N, VD(N)) ;    // triangular superior
    VD P(N);               // vector de permutacion
    MD A_Copy(N, VD(N));   // copia de la matriz
    A_Copy = A;
    
    // vector de perm dim N+1, posicion N para numero de operaciones
    for (i = 0; i <= N; i++) P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        int imax = find_max_pivot(A, N, i); // pos pivote max, parcial escalonado
        if (abs(A[imax][i]) < tol) return 0; //failure, matrix is degenerate
        swap(A[i], A[imax]);
        swap(P[i], P[imax]);
        P[N]++;

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i]; // Triangular inferior

            for (k = i + 1; k < N; k++){
                A[j][k] -= A[j][i] * A[i][k];
            }
        }
    }
    cout << "this is A, contains L and U" << endl;
    write(A, N);
    cout << endl;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j) U[i][i] = A[i][i];
            else if (i > j) L[i][j] = A[i][j];
            else U[i][j] = A[i][j];
        }
    }
    
    // ** Calculo del determinante ** //
    double det = 1.;
    for (int i = 0; i < N; i++) det *= U[i][i];
    
    cout << "this is U" << endl;
    write(U, N);
    cout << "determinante: " << det << endl; 
    
    LUPSolve(A, P, b, N);
    
    return 1;  //decomposition done
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
void LUPSolve(MD& A, VD& P, VD& b, int N) {

    VD x(N);
    
    for (int i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (int k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];

        x[i] = x[i] / A[i][i];
    }
    
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


void write(MD& A, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl; cout << endl; cout << endl;
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


