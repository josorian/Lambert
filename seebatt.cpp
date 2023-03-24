//
// Created by josorian on 22/03/2023.
//

#include "seebatt.h"

double seebatt(double v) {
    double c[21] = {0.2, 9.0 / 35.0, 16.0 / 63.0, 25.0 / 99.0, 36.0 / 143.0, 49.0 / 195.0, 64.0 / 255.0, 81.0 / 323.0,
                    100.0 / 399.0, 121.0 / 483.0, 144.0 / 575.0, 169.0 / 675.0, 196.0 / 783.0, 225.0 / 899.0,
                    256.0 / 1023.0, 289.0 / 1155.0, 324.0 / 1295.0, 361.0 / 1443.0, 400.0 / 1599.0, 441.0 / 1763.0,
                    484.0 / 1935.0};
    double sqrtopv = sqrt(1.0 + v);
    double eta = v / pow(1.0 + sqrtopv, 2);
    double delold = 1.0;
    double temold = c[0];
    double sum1 = temold;
    double del;
    double term;
    int i = 0;
    while (i <= 20 && fabs(temold) > 0.00000001) {
        del = 1.0 / (1.0 + c[i + 1] * eta * delold);
        term = temold * (del - 1.0);
        sum1 = sum1 + term;
        i++;
        delold = del;
        temold = term;
    }
    return 1.0 / ((1.0 / (8.0 * (1.0 + sqrtopv))) * (3.0 + sum1 / (1.0 + eta * sum1)));
}
