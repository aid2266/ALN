#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
//#include "lu.h"
//#include "resol.h"
//#include <fstream> // utilizamos para leer files

using namespace std;

const double tol = 0.0000000000001;

typedef vector<double> VD;
typedef vector <VD> MD;

void lu(MD& A, int n, double tol);
void write (MD& A, int n);

int main(){
    
    int n;
    cin >> n; // size of matrix
    MD A(n, VD(n)); // init A
    VD b(n); // init b
    
    cerr << "Enter matrix A of size " << n << "x" << n << endl;
    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) cin >> A[i][j]; // reading the matrix
    }
    
    cerr << "Enter matrix b " << endl;
    for (int i = 0; i < n; i++) cin >> b[i]; // reading b
    
    lu(A, n, tol); // call the function
    
}

static bool abs_compare(double a, double b){
    return (abs(double(a)) < abs(double(b)));
}


void lu(MD& A, int n, double tol){
    
    MD L(n, VD(n)); // lower triangular
    MD U(n, VD(n)) ; // upper triangular
    VD perm(n); // vector de perm
    
    cerr << "you will LU factorise the following matrix: " << endl;
    write(A, n);
    
    for (int i = 0; i < n; i++) perm[i] = i; // init vector perm
    
    for (int k = 0; k < n; k++){
        //U[k][k] = A[k][k]; // choosing pivot
        L[k][k] = 1; // set diagonals to 1
        VD s(n); // creamos un vector de tamaño
        
        /*
        for (int i = k; i < n; i++){
            s[i] = *max_element(A[i].begin(), A[i].end(), abs_compare); // elemento max
            double temp = distance(A[i].begin(), A[i].end());
            cout << "max element at: " << temp << endl; 
        }
        
        double temp = 0;
        for (int i = 0; i < n; i++){
            if (temp < A[i][k]/s[i]) temp  = A[i][k]/s[i];
        }
        
        cout << temp << endl;
        */
        
        for (int i = k+1; i < n; i++){
            L[i][k] = A[i][k] / U[k][k];
            U[k][i] = A[k][i];
        }
        for (int i = k+1; i < n; i++){
            for (int j = k+1 ; j < n; j++){
                A[i][j] = A[i][j] - L[i][k]*U[k][j];
            }
        }
    }
    write(L, n);
    write(U, n);
}


void write(MD& A, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}
