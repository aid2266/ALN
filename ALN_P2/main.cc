#include <iostream>
#include <vector>
#include "lu.h"
#include "resol.h"

//#include <fstream> // utilizamos para leer files
const double tol = 0.0000000000001;


using namespace std;

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
    
    // llegeix la linea de comandes -- fin , fout
    // fin <- n
    // reserva dinamica de memoria, para matrices A y b
    // FEM copia de matriu A per calcular el error
    // utilizar pointers o STL
    
    // CRIDEM LU per fer el calcul del error
    // necessitem la matriu original (abans de fer LU)
    // ABANS de LU fem copia de A, per veure error
    // CALCULEM ERROR per norma 1 de PA - LU
    
    // CALCULEM ERROR residu - norma 1 de |Bx - b| (on B es copia de A, x sol)

    // OUT PER LA PANTALLA: 
    // - size n de la matriu
    // - residu de la descomposicio  (error)
    // - vector perm
    // - determinant de la matriu
    // - residuu de solucio Bx - b (per les tres normes)
    
    

}


// DESCOMPOSICIO LU
vector< vector<double> > lu(MD& A, vector<double>& b, double tol){
    int n = int(A.size()) // size of matrix A
    // RECORDAMOS: queremos hacer LU encima de la matriz
    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            double pivot = A[i][i];
            double m = A[i][j] / pivot; // multiplicador
            A[i][j] -= m*A;
        }
    }
    
    
}
