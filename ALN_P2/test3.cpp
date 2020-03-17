#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

typedef vector<double> VD;
typedef vector <VD> MD;

void write(MD& A, int n);

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
    
    for (int counter = 0; counter < v; counter++){
        int i;
        double elem;
        inFile >> i >> elem;
        b[i] = elem;
    }
    
    write(A, n);
    
    return 0;
    
}


void write(MD& A, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A[i][j] << ' ';
        }
        cout << endl;
    }
}
