#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// empleamos el algoritmo de horner

int main(){
    
    double n;
    double x;
    cin >> n >> x; 
    double p = 1; // init a 1 
    
    for (double i = n; i > 0; i--){
        p = p * x/i + 1; 
    }
    
    cout << p << endl; 
    
    /*
    cerr << "enter N: " << endl;
    cin >> n;
    cerr << "enter value of x: " << endl;
    cin >> x;
    cout << Maclaurin(n,x) << endl;
    //cout << setprecision(4) << sum << endl;
    * */
    /*
    //vector<double> VD(n); //longitd n
    //float p = VD[n]; // init a la pos n-1
    for (int i = n; i > 0; i--){
        p = p * x + a[i-1]; 
    }
    */ 
}
