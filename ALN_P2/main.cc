#include <iostream>
#include <vector>
#include "lu.h"
#include "resol.h"
//#include <fstream> // utilizamos para leer files

using namespace std;

const double tol = 0.0000000000001;

typedef vector<double> VD;
typedef vector <VD> MD;


int main(){
    
    int n;
    cin >> n; // size of matrix
    vector< vector<double> > A(n, vector<double>(n)); // init matrix A of size nxn
    vector<double> b(n); // init vector solution b
    
    cerr << "Enter matrix A of size " << n << "x" << n << endl;
    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) cin >> A[i][j]; // reading the matrix
    }
    
    cerr << "Enter matrix b " << endl;
    for (int i = 0; i < n; i++) cin >> b[i]; // reading b
    
    // vector< vector<double>> cA = A; // make a copy of matrix A
    
    lu(A, b, tol);

}

/*
vector< vector<double> > lu(MD& A, vector<double>& b, double tol){
    int n = int(A.size()) // size of matrix A
    // RECORDAMOS: idealmente queremos hacer LU encima de la matriz
    // DECLARAMOS: dos matrices, L y U
    MD(n, VD(n, 0)) L; // lower triangular
    MD(n, VD(n, 0)) U; // upper triangular
    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            double pivot = A[i][i];
            double m = A[i][j] / pivot; // multiplicador
            A[i][j] -= m*A;
        }
    }
    
    
}
*/
