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


double find_max_pivot (const MD& A, int n, int k){
    double pos_max = k; // index of the row with max pivot
    double max_pivot = abs(A[k][k]); // assume it is initial pivot
    
    for (int i = k; i < n; i++){
        
        double elem = abs(A[i][k]); // current element of column
        if (elem > max_pivot){
            max_pivot = elem;
            pos_max = i;
        }
    }
    
    cerr << "position max is: " << pos_max << endl;
    
    return pos_max;
    
}

void lu (MD& A, int n, double tol){
    
    VD perm(n);
    for (int i = 0; i < n; i++) perm[i] = i; // init perm matrix
    
    for (int k = 0; k < n; k++){
        double pos_max = find_max_pivot(A, n, k);
        if (abs(A[pos_max][k]) < tol) cout << "la matriz es singular";
        
        swap(A[k], A[pos_max]); // hacer un swap
        swap(perm[k], perm[pos_max]);
        
    }
    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A[i][j] << ' ';
        }
        cout << endl;
    }
    
}
