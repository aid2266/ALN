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

void lu(MD& A, int n, double tol);
void write (MD& A, int n);
void writeV (VD& b, int n);
double find_max_pivot (const MD& A, int n, int k);
static bool abs_compare(double a, double b);

int main(){
    
    ifstream inFile;
    inFile.open("test.txt");
    
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
    
    lu(A, n, tol); // call the function
    
    // ** una vez realizada la descomposicion LU ** //
    cout << "dimension del sistema: " << n << endl;
    cout << "estimacion del error de descomposicion PA = LU: " << 'e' << endl;
    cout << "vector de permutaciones de P: " << 'p' << endl;
    cout << "determinant de la matriu A: " << 'd' << endl;
    cout << "estimacio del error del sistema Ax = b, amb norma infinita: " << 'i' << endl;
    
}


// *** function that does LU factorisation ** //

void lu(MD& A, int n, double tol){
    
    MD L(n, VD(n)); // lower triangular
    MD U(n, VD(n)) ; // upper triangular
    VD perm(n); // vector de perm
    
    for (int i = 0; i < n; i++) perm[i] = i; // init vector perm
    
    // ** LU FACTORISATION ** //
    for (int k = 0; k < n; k++){
        L[k][k] = 1; // set diagonal of L
        
        double pos_max = find_max_pivot(A, n, k);
        if (abs(A[pos_max][k]) < tol) cout << "la matriz es singular"; // tecnicamente devuelve -1
        
        swap(A[k], A[pos_max]); // hacer un swap
        swap(perm[k], perm[pos_max]);
        
        //cerr << "matrix A after permutation: " << endl;
        //write(A, n);
        //cerr << endl;
        
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
        
        //cerr << "this is the k:" << k << " iteration, with matrix A:" << endl;
        //write(A, n);
        
    }
        
    cerr << "print vector de perm: " << endl;
    for (int i = 0; i < n; i++) cout << perm[i] << ' ';
    cout << endl;
    cerr << "matrix A: " << endl;
    write(A, n);
    cerr << "matrix L: " << endl;
    write(L, n);
    cerr << "matrix U: " << endl;
    write(U, n);
    cerr << "matrix perm: " << endl;
    writeV(perm, n);

    
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
    for (int j = 0; j < n; j++){
            cout << b[j] << ' ';
        }
    cout << endl;
}

// *** function that finds max pivot *** //

double find_max_pivot (const MD& A, int n, int k){
    double pos_max = k; // index of the row with max pivot
    double max_pivot = 0; // assume it is initial pivot
    
    for (int i = k; i < n; i++){
        
        double temp = *max_element(A[i].begin(), A[i].end(), abs_compare); // max element of row
        //cerr << "the maximum value of the row is: " << temp << endl;
        
        double elem = abs(A[i][k]/temp); // current element of column
        cerr << "current element is " << elem << endl;
        if (elem > max_pivot){
            max_pivot = elem;
            pos_max = i;
        }
    }
    
    cerr << "position max is: " << pos_max << endl;
    
    return pos_max;
    
}

// ** Function that compares absolute values results ** //
static bool abs_compare(double a, double b){
    return (abs(double(a)) < abs(double(b)));
}
