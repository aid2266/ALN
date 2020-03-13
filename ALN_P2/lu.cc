#include "lu.h" // where we include the file
//#include <iostream>
//#include <vector>

using namespace std; 

// LU DECOMPOSITION FILE

void lu(){
    cout << "Hello world" << endl;
    return;
}
    // tolerancia - 10^-12
    // utilitzem la tolerancia per veure si la matriu es singular
    // creem una funcio que fa la descomposicio LU
    // la funcio retorna en LA MATEIXA matriu la descomp LU
    // tambe retorna la matriu permutacio per referencia

    // si hem fet pivotatge, tenim vector perm
    // retornem 1 si descomposicio OK, nombre parell de perm
    // retornem -1 si descomposicio OK, pero nombre senar perm
    // retornem 0 si matriu es singular
    
    // CALCUL DE ERROR: norma 1 de PA - LU (NO!)
    // fem el calcul a main

