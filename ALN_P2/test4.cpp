#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
//#include "lu.h"
//#include "resol.h"

using namespace std;

const double tol = 0.0000000000001;

typedef vector<double> VD;
typedef vector <VD> MD;

int lu(MD& A, int n, double tol);
void write (MD& A, int n);
void writeV (VD& b, int n);
MD multiplyMatrix(MD& A, MD& B, int n);
double find_max_pivot (const MD& A, int n, int k);
static bool abs_compare(double a, double b);

int main(){
    
    ifstream inFile;
    inFile.open("/Users/aidandeaves/Documents/ALN/ALN_P1/ALN_P2/MAT/M01.DAT"); // change this file!!
    
    // check for any errors
    if (inFile.fail()){
        cerr << "Error Opening File" << endl;
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
    
    cout << "your matrix A is the following: " << endl;
    write(A, n);
    cout << "your matrix b is the following: " << endl;
    writeV(b, n);
    
    cout << lu(A, n, tol) << endl; // call the function
    
    // ** una vez realizada la descomposicion LU ** //
    cout << "dimension del sistema: " << n << endl;
    cout << "estimacion del error de descomposicion PA = LU: " << 'e' << endl;
    cout << "vector de permutaciones de P: " << 'p' << endl;
    cout << "determinant de la matriu A: " << 'd' << endl;
    cout << "estimacio del error del sistema Ax = b, amb norma infinita: " << 'i' << endl;
    
}


// *** function that does LU factorisation ** //

int lu(MD& A, int n, double tol){
    
    MD L(n, VD(n));     // lower triangular
    MD U(n, VD(n)) ;    // upper triangular
    MD A_Copy(n, VD(n));// copy of matrix A
    VD perm(n); // vector de perm
    int numPermutations = 0;
    
    for (int i = 0; i < n; i++) perm[i] = i; // init vector perm
    A_Copy = A;         // make a copy of A
    
    // ** LU FACTORISATION ** //
    for (int k = 0; k < n; k++){
        L[k][k] = 1; // set diagonal of L
        double pos_max = find_max_pivot(A, n, k);
        if (abs(A[pos_max][k]) < tol) return 0;
        
        swap(A[k], A[pos_max]); // hacer un swap
        swap(perm[k], perm[pos_max]);
        numPermutations++;
        
        U[k][k] = A[k][k]; // este valor es el que cogemos como pivote!!
        
        for (int i = k+1; i < n; i++){
            L[i][k] = A[i][k] / U[k][k];
            U[k][i] = A[k][i];
        }
        for (int i = k+1; i < n; i++){
            for (int j = k+1 ; j < n; j++){
                A[i][j] = A[i][j] - L[i][k]*U[k][j]; // renombramos la matriz A
            }
        }
    }
        
    // ** CALCUL ERROR PA = LU ** //
    // ** CALCUL DE LU ** //
    MD C(n, VD(n));
    C = multiplyMatrix(L, U, n);
    
    // ** CALCUL DEL ERROR ** //
    // double norma_1 = 0;
    double norma_1 = 0;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            /*
            double LU_ij = 0; // valor de LU en [i][j]
            for (int k = 0; k < n; k++){
                LU_ij += L[i][k] * U[k][j]; // esto calcula L*U en casilla [i][j]
            }
            norma_1 += A_Copy[perm[i]][j] - LU_ij;
            */
            norma_1 += A_Copy[perm[i]][j] - C[i][j];
        }
    }
    
    // ** CALCULO DEL DETERMINANTE ** //
    double det = 1;
    for (int i = 0; i < n; i++){
        det *= U[i][i];
    }
    
    
    cerr << "matrix A: " << endl;
    write(A_Copy, n);
    cerr << "matrix L: " << endl;
    write(L, n);
    cerr << "matrix U: " << endl;
    write(U, n);
    cerr << "this is the matrix LU" << endl;
    write(C, n);
    cout << "el error |PA - LU| es: " << abs(norma_1) << endl;
    cout << "el vector permutación es: ";
    writeV(perm, n);
    cout << "check # permutaciones, si matriz es singular" << endl;
    cout << "detLU: " << det << endl;
    cerr << "num of permutations " << numPermutations << endl;
    
    // ** ||PA - LU||
    if (numPermutations % 2 == 0) return 1;
    else return -1;
    
}

// ** function that multiplies two matrices ** //
MD multiplyMatrix(MD& A, MD& B, int n){
    MD C(n, VD(n)); // empty matrix
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            for (int k = 0; k < n; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}


// *** function that writes matrix *** //
void write(MD& A, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}
// *** function that writes vector *** //
void writeV(VD& b, int n){
    cout << '[';
    for (int j = 0; j < n; j++){
            cout << b[j] << ' ';
        }
    cout << ']' << endl;
}

// *** function that finds max pivot *** //

double find_max_pivot (const MD& A, int n, int k){
    double pos_max = k; // index of the row with max pivot
    double max_pivot = 0; // assume it is initial pivot
    
    for (int i = k; i < n; i++){
        
        double temp = *max_element(A[i].begin(), A[i].end(), abs_compare); // max element of row
        //cerr << "the maximum value of the row is: " << temp << endl;
        
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
