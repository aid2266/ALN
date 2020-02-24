#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double Maclaurin (double k, double x){
    if (k == 0) return 1;
    //else return (1/k)*x + Maclaurin(k-1, x);
    else return Maclaurin (k-1, x)/k + x*Maclaurin(k-1, x)/k;
}


int main(){
    
    double n;
    double x;
    cerr << "enter N: " << endl;
    cin >> n;
    cerr << "enter value of x: " << endl;
    cin >> x;
    cout << Maclaurin(n,x) << endl;
    //cout << setprecision(4) << sum << endl;
}
