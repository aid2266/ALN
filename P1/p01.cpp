#include <iostream>
#include <cmath>
#include <cfloat>


// CAMBIOS: set float epsilon a 1
// importar cfloat para soluciones ya hechas

using namespace std;

int main(){
    
    float ep1;
    double ep2;
    cin >> ep1;
    ep2 = ep1;
    while (float(1 + ep1) > 1) ep1 /= 2;
    while (double(1 + ep2) > 1) ep2 /= 2;
    cout << "Float value of epsilon: " << 2 * ep1 << endl;
    cout << "Double value of epsilon: " << 2 * ep2 << endl;
    
    cout << "Real values double: " << FLT_EPSILON << endl; 
    cout << "Real values double: " << DBL_EPSILON << endl; 
    
    // comandos con cfloat 
    // FLT_EPSILON
    // DBL_EPSILON
}
