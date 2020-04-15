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

typedef vector<double> VD;  // vector de doubles
typedef vector <VD> MD;     // matriz de doubles

// creamos headers que son llamados por main //
int lu(MD& A, VD& b, VD& x, int n, double tol);
void resol(MD& A, MD& A_Copy, VD& P, VD& b, VD& x, int N);


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
    VD x(n);            // vector solucion x
    
    for (int counter = 0; counter < v; counter++){
        int i, j;
        double elem;
        inFile >> i >> j >> elem;
        A[i][j] = elem; // leemos la matriz del fichero
    }
    
    inFile >> k;
    
    for (int counter = 0; counter < k; counter++){
        int i;
        double elem;
        inFile >> i >> elem;
        b[i] = elem;
    }
    
    inFile.close(); // cerramos inFile para prevenir errores
    cout << lu(A, b, x, n, tol) << endl; // llamamos funcion lu
    
    ofstream outFile; // recordamos que tenemos x pasado por referencia
    outFile.open("output.txt"); // creamos fichero (o rewrite si ya esta creado)
    for (int i = 0; i < n; i++){
        outFile << i << " " << x[i] << endl; //
    }
    outFile.close();
    
}
