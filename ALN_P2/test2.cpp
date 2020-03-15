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
void write(MD& A, int n);

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


double find_max_pivot (const MD& A, int n, int k){
    double pos_max = k; // index of the row with max pivot
    double max_pivot = 0; // assume it is initial pivot
    
    for (int i = k; i < n; i++){
        
        double temp = *max_element(A[i].begin(), A[i].end()); // max element of row
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

void lu (MD& A, int n, double tol){
    
    VD perm(n); // init vector de permutacion
    for (int i = 0; i < n; i++) perm[i] = i; // init perm matrix
    
    cerr << "you are going to diagonalise the following matrix A: " << endl;
    write(A, n);
    
    MD L(n, VD(n)); // lower triangular
    MD U(n, VD(n)) ; // upper triangular
    
    for (int k = 0; k < n; k++){
        L[k][k] = 1; // set diagonal of L
        
        double pos_max = find_max_pivot(A, n, k);
        if (abs(A[pos_max][k]) < tol) cout << "la matriz es singular";
        
        swap(A[k], A[pos_max]); // hacer un swap
        swap(perm[k], perm[pos_max]);
        
        U[k][k] = A[k][k]; // este valor es correcto!!
        
        for (int i = k+1; i < n; i++){
            L[i][k] = A[i][k] / U[k][k];
            U[k][i] = A[k][i];
        }
        for (int i = k+1; i < n; i++){
            for (int j = k+1 ; j < n; j++){
                A[i][j] = A[i][j] - L[i][k]*U[k][j];
            }
        }
        
        cerr << "this is the k:" << k << " iteration, with matrix A:" << endl;
        write(A, n); 
        
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
