#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// P003 NUMERO AUREO DE TRES METODOS

const float phi = 1.618033988749;

float phi1 (int k){
    if (k == 0) return 1;
    if (k == 1) return phi;
    else return phi * phi1(k-1);
}

float phi2 (int k){
    if (k == 0) return 1;
    if (k == 1) return phi;
    if (k == 2) return float (1 - phi);
    else return phi2(k - 2)*(float(1 - phi));
}

float phi3 (int k){
    if (k == 0) return 1;
    if (k == 1) return phi;
    if (k == 2) return float (1- phi);
    else return phi3(k - 2) - phi3(k - 1);
}

int main(){
    
    int num;
    cin >> num; // this is the constant of power
    cout << "Phi Method 1: " << phi1(num) << endl;
    cout << "Phi Method 2: " << phi2(num) << endl;
    cout << "Phi Method 3: " << phi3(num) << endl;
}
