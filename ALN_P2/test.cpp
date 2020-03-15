#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

const double tol = 0.0000000000001;
const double epsilon = 0.0000000000001;
typedef vector<double> VD;
typedef vector <VD> MD;

void lu(MD& A, int n, double tol);

int main(){
    
    int n;
    cin >> n; // size of matrix
    MD A(n, VD(n)); // init A
    cerr << "Enter matrix A of size " << n << "x" << n << endl;
    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) cin >> A[i][j]; // reading the matrix
    }
    
    lu(A, n, tol); // call the function
    
}

static bool abs_compare(double a, double b){
    return (abs(double(a)) < abs(double(b)));
}


void lu (MD& A, int n, double tol){
    
    // queremos el maximo de cada fila
    for (int i = 0; i < n; i++){
        MD s(n);
        for (int j = 0; j < n; j++){
            s[j] = *max_element(A[j].begin(), A[j].end(), abs_compare);
        }
        // ahora tenemos el elemento maximo de cada fila!
        double pivot_max = abs(A[i][i]);
        double pos_max = 0;
        
        for (int k = i; k < n; k++){
            
            cerr << "en el paso k = " << k << " tenemos: " << A[k][i] << '/' << s[k] << endl;
            double temp = abs(A[k][i]/s[k]);
            if (temp > pivot_max) {
                pivot_max = temp;
                pos_max = k; // k representa la distancia
            }
        }
        
        // ahora tenemos pivot max y pos max
        cerr << "el calculo maximo del pivote es: " << pivot_max << endl;
        cerr << "el calculo de la posicion de este pivot es: " << pos_max << endl;
        
    }
    
}
