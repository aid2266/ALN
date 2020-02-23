#include <iostream>
#include <cmath>

using namespace std;

int main(){
    
    float ep1;
    double ep2;
    cin >> ep1;
    ep2 = ep1;
    while (float(1 + ep1) > 1) ep1 /= 2;
    while (double(1 + ep2) > 1) ep2 /= 2;
    cout << "Float value of epsilon: " << ep1 << endl;
    cout << "Double value of epsilon: " << ep2 << endl;

}
