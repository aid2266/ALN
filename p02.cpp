#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

int main(){
    
    float sumUp = 0;
    float sumDown = 0;
    for (float i = 1; i <= 15; i++) sumUp += 1 / (i * i);
    cout << "Sum Up: " << setprecision(4) << sumUp << endl;
    // como escoger por eliminacion o sin eliminacion?
    
    for (float i = 15; i >= 1; i--) sumDown += 1 / (i * i);
    cout << "Sum Down: " << setprecision(4) << sumDown << endl;

}
