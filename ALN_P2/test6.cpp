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

typedef vector<double> VD;
typedef vector <VD> MD;
int lu(MD& A, VD& b, int n, double tol);
int LUPDecompose(MD& A, VD& b, int n, double tol);
int copyLUP(MD& A, VD& b, int n, double tol);
void resol(MD& L, MD& U, VD& b, VD& perm, MD& A, int n);
MD multiplyMatrix(MD& A, MD& B, int n);
int find_max_pivot (const MD& A, int n, int k);
static bool abs_compare(double a, double b);
void write (MD& A, int n);
void writeV (VD& b, int n);

int main(){
    
    ifstream inFile;
    inFile.open("/Users/aidandeaves/Documents/ALN/ALN_P1/ALN_P2/MAT/M02.DAT"); // change this file!!
    
    // check for any errors
    if (inFile.fail()){
        cerr << "Error Obrint Fitxer" << endl;
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
    
    /*
    cout << "your matrix A is the following: " << endl;
    write(A, n);
    cout << "your matrix b is the following: " << endl;
    writeV(b, n);
    */
    
    cout << copyLUP(A, b, n, tol) << endl; // call the function
    
    // ** una vez realizada la descomposicion LU ** //
    cout << "dimension del sistema: " << n << endl;
    cout << "estimacion del error de descomposicion PA = LU: " << 'e' << endl;
    cout << "vector de permutaciones de P: " << 'p' << endl;
    cout << "determinant de la matriu A: " << 'd' << endl;
    cout << "estimacio del error del sistema Ax = b, amb norma infinita: " << 'i' << endl;
    
}

int copyLUP(MD& A, VD& b, int n, double tol){
    
    MD L(n, VD(n));     // triangular inferior
    MD U(n, VD(n)) ;    // triangular superior
    MD A_Copy(n, VD(n));// copy of matrix A
    VD perm(n);         // vector de perm
    int numPermutations = 0;
    
    for (int i = 0; i < n; i++) perm[i] = i; // init vector perm
    A_Copy = A;         // make a copy of A
    
    cerr << "this is the matrix you will LU decompose: " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A_Copy[i][j] << ' ';
        }
        cout << '|' << b[i];
        cout << endl;
    }
    
    
    // ** LU FACTORISATION ** //
    for (int k = 0; k < n; k++){
        L[k][k] = 1; // set diagonal of L
        int pos_max = find_max_pivot(A, n, k); // pos pivote max, parcial escalonado
        //if (abs(A[pos_max][k]) < tol) return 0; // matriz es singular
        
        swap(A[k], A[pos_max]); // permutamos filas
        swap(b[k], b[pos_max]); // permutamos b
        swap(perm[k], perm[pos_max]); // modificamos matriz perm
        numPermutations++; // actualizamos numero de permutaciones
        
        U[k][k] = A[k][k];
        
        // triangular superior
        for (int j = k; j < n; j++){
            double sum = 0.;
            for (int p = 0; p < k; p++) sum += L[k][p]*U[p][j];
            U[k][j] = A[k][j] - sum; // triangular inferior
        }
        // triangular inferior
        for (int i = k+1; i < n; i++){
            double sum = 0.;
            for (int p = 0; p < k; p++) sum += L[i][p]*U[p][k];
            L[i][k] = (A[i][k] - sum) / U[k][k];
        }
    }
    
    MD C(n, VD(n));
    C = multiplyMatrix(L, U, n);
    
    cerr << "This is the matrix L" << endl;
    write(L, n);
    cout << endl;
    cout << endl;
    cerr << "this is the matrix U" << endl;
    write(U, n);
    cout << endl;
    cout << endl;
    cerr << "this is the matrix LU with permutations" << endl;
    write(C, n);
    cout << endl;
    cout << endl;
    cerr << "this is matrix PA" << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A_Copy[perm[i]][j] << ' ';
        }
        cout << endl;
    }
    cout << endl; cout << endl; 
        
    // ** CALCUL ERROR PA = LU ** //
    // ** CALCUL DE LU ** //
    
    // ** CALCUL DEL ERROR ||PA = LU||_1 ** //
    double norma_1 = 0;
    //double test_norma = 0.;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            norma_1 += A_Copy[perm[i]][j] - C[i][j]; // Copia de A permutada
            //test_norma += A_Copy[i][j] - C[i][j];
        }
    }
    
    // ** CALCULO DEL DETERMINANTE ** //
    double det = 1.;
    for (int i = 0; i < n; i++){
        det *= U[i][i];
    }

    cerr << "matrix A after permutations with b: " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A_Copy[perm[i]][j] << ' ';
        }
        cout << '|' << b[i];
        cout << endl;
    }
    
    cerr << "this is permuted b " << endl;
    for (int i = 0; i < n; i++) cout << b[perm[i]] << ' ';
    cout << endl;
    cerr << "this is also permuted b" << endl;
    writeV(b, n);
    
    cout << "el error |PA - LU| es: " << setprecision(10) << abs(norma_1) << endl;
    //cout << "el error test |PA - LU| es" << abs(test_norma) << endl;
    cout << "el vector permutación es: ";
    writeV(perm, n);
    cout << "detLU: " << det << endl;
    cerr << "num of permutations " << numPermutations << endl;
    
    resol(L, U, b, perm, A_Copy, n); // resuelve la matriz
    
    // ** ||PA - LU||
    if (numPermutations % 2 == 0) return 1;
    else return -1;
    
}



// *** function that does LU factorisation ** //

void resol(MD& L, MD& U, VD& b, VD& perm, MD& A, int n){
    VD y(n); // creamos vector y aux
    VD x(n); // vector solucion
    
    // ** resolvemos Ly = b ** //
    for (int i = 0; i < n; i++){
        y[i] = b[perm[i]];
        for (int j = 0; j < i; j++){
            y[i] = y[i] - L[i][j]*y[j];
        }
    }
    cerr << "this is vector y: " << endl;
    writeV(y, n);
    
    // ** resolvemos Ux = y ** //
    for (int i = n-1; i >= 0; i--){
        x[i] = y[i];
        for (int j = i+1; j < n; j++){
            x[i] = x[i] - U[i][j]*x[j];
        }
        x[i] = x[i] / U[i][i];
    }
    /*
    for (int i = n-1; i >= 0; i--){
        double sum = 0;
        for (int j = n-1; j > i; j--) sum += U[i][j] * x[j];
        x[i] = (y[i] - sum) / U[i][i];
    }
    */
    cerr << "this is solution vector x: " << endl;
    writeV(x, n);
    
    // CALCULO DE LA NORMA ||Ax - b|| aprox 0
    double norma_1 = 0.;
    double norma_2 = 0.;
    double norma_infinity = 0;
    for (int i = 0; i < n; i++){
        double sum_1 = 0.;
        double sum_infinity = 0;
        for (int j = 0; j < n; j++){
            sum_1 += A[i][j]*x[j]; // valor de b'
        }
        norma_1 += sum_1 - b[i]; //
        norma_2 += (sum_1 - b[i])*(sum_1 - b[i]);
    }
    
    cout << "este es el ERROR Ax - b norma 1 : " << norma_1 << endl;
    cout << "este es el ERROR Ax - b norma 2 : " << norma_2  << endl;

}

// BORRAR //

int LUPDecompose(MD& A, VD& b, int n, double tol) {


    MD A_Copy;
    VD perm(n);             // vector de perm
    int numPermutations = 0;
    
    for (int i = 0; i < n; i++) perm[i] = i; // init vector perm
    A_Copy = A;         // make a copy of A

    for (int i = 0; i < n; i++) {
        int pos_max = find_max_pivot(A, n, i); // pos pivote max, parcial escalonado
        if (abs(A[pos_max][i]) < tol) return 0; // matriz es singular
        
        swap(A[i], A[pos_max]); // permutamos filas
        //swap(b[k], b[pos_max]); // permutamos b
        swap(perm[i], perm[pos_max]); // modificamos matriz perm
        numPermutations++; // actualizamos numero de permutaciones

        for (int j = i + 1; j < n; j++) {
            A[j][i] /= A[i][i]; // multiplicador

            for (int k = i + 1; k < n; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }
    cerr << "this is your decomposed LUP matrix" << endl;
    write(A, n);

    return 1;  //decomposition done
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
    cout << endl; cout << endl; cout << endl;
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

int find_max_pivot (const MD& A, int n, int k){
    int pos_max = k; // index of the row with max pivot
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

