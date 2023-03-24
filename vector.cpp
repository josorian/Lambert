//
// Created by josorian on 15/03/2023.
//

#include "vector.h"

//Documentar
double norm(vector<double> v, int n) {
    double suma = 0;
    if (n <= 0) {
        throw "Empty vector";
    }
    for (int i = 0; i < n; i++) {
        suma += v[i] * v[i];
    }
    return sqrt(suma);
}

double dot(vector<double> v1, vector<double> v2, int n1, int n2) {
    double suma = 0;
    if (n1!=n2) {
        throw "size diferent";
    }
    for(int i=0;i<n1;i++){
        suma +=v1[i]*v2[i];
    }
    return suma;
}

vector<double> cross(vector<double> v1, vector<double>v2, int n1, int n2) {
    if (n1!=n2) {
        throw "size diferent";
    }
    if (n1 != 3 & n2 != 3){
        throw "n1 or n2 diferent 3";
    }
    vector<double> resul(3);
    resul[0]=v1[0]*v2[0]-v1[2]*v2[1];
    resul[1]=v1[2]*v2[0]-v1[0]*v2[2];
    resul[2]=v1[0]*v2[1]-v1[1]*v2[0];
    return resul;
}

