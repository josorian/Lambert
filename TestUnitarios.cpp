#include <iostream>
#include <cmath>
#include "vector.h"
#include "seebatt.h"
#include "seebattk.h"

using namespace std;

int main() {
    vector<double> v = {1, 1, 1};
    if (norm(v, 3) - sqrt(3.0) < pow(10, -12)) {
        cout << "Test norm(): passed" << endl;
    } else {
        cout << "Test norm(): failed" << endl;
    }
    vector<double> v1 = {1, 1, 2};
    vector<double> v2 = {1, 2, 3};
    if (dot(v1, v2, 3, 3) - (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) < pow(10, -12)) {
        cout << "Test dot(): passed" << endl;
    } else {
        cout << "Test norm(): failed" << endl;
    }
    vector<double> c = cross(v1, v2,3,3);
    if (c[0] - (v1[0] * v2[0] - v1[2] * v2[1]) < pow(10, 12) &&
        c[1] - (v1[2] * v2[0] - v1[0] * v2[2]) < pow(10, 12) & c[2] - (v1[0] * v2[1] - v1[1] * v2[0]) < pow(10, 12)) {
        cout << "Test cross(): passed" << endl;
    } else {
        cout << "Test norm(): failed" << endl;
    }
    if (seebatt(0.0) - 5.0 <= pow(10, -14)) {
        cout << "seebatt(0.0) es correcto" << endl;
    } else {
        cout << "Falla: seebatt(0.0)" << endl;
    }
    if (seebatt(1.0) - 6.06251 <= pow(10, -5)) {
        cout << "seebatt(0.0) es correcto" << endl;
    } else {
        cout << "seebatt(1.0) falla" << endl;
    }
    if (seebatt(10000.0) - 257.145 <= pow(10, -10)) {
        cout << "seebatt(0.0) es correcto" << endl;
    } else {
        cout << seebatt(10000.0) << endl;
    }
    if (seebattk(0.0) - 0.333333 <= pow(10,-5)){
        cout<<"seebattk(0.0) es correcto"<<endl;
    }else{
        cout<<"seebattk(0.0) falla"<<endl;
    }
    if (seebattk(0.0) - 0.333333 <= pow(10,-5)){
        cout<<"seebattk(0.0) es correcto"<<endl;
    }else{
        cout<<"seebattk(0.0) falla"<<endl;
    }
    if (seebattk(100.0) - 0.0211 <= pow(10,-5)){
        cout<<"seebattk(100.0) es correcto"<<endl;
    }else{
        cout<<"seebattk(100.0) falla"<<endl;
    }
}