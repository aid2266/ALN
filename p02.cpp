#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// SUMA FINITA DE 15 TERMINOS DE 1/n^2
// queremos encontrar el error en los calculos!! 

// RESULTADOS:  0.153 x 10 ^ 1
//              0.157 x 10 ^ 1

// COMO TRUNCAR EL ERROR? 
// exemple - X = 10.33251, ho vull troncar en 5 xsf, t = 5
//          passem a exponencial : 0.1044251 X 10 ^ 2 
//          int(log10 abs(x)) + 1, esto es el exponente = 2
//          x * 10 ^(-2 + t) on t es el nombre de xsf
//          10332.51 vull troncar,C++ : trunc
//          trunc(x) -- 10322.51 x 10 ^ (2-5) = 0.10322 x 10 ^ 2

//          si volem per aproximacio -- round trunc  


int main(){
    
    float sumUp = 0;
    float sumDown = 0;
    for (float i = 1; i <= 15; i++) sumUp += 1 / (i * i);
    cout << "Sum Up: " << setprecision(4) << sumUp << endl;
    // como escoger por eliminacion o sin eliminacion?
    
    for (float i = 15; i >= 1; i--) sumDown += 1 / (i * i);
    cout << "Sum Down: " << setprecision(4) << sumDown << endl;

}
