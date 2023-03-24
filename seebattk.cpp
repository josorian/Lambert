//
// Created by josorian on 22/03/2023.
//

#include "seebattk.h"

double seebattk(double v) {
    double d[21] = {1.0 / 3.0, 4.0 / 27.0, 8.0 / 27.0, 2.0 / 9.0, 22.0 / 81.0, 208.0 / 891.0,
                    340.0 / 1287.0, 418.0 / 1755.0, 598.0 / 2295.0, 700.0 / 2907.0, 928.0 / 3591.0,
                    1054.0 / 4347.0, 1330.0 / 5175.0, 1480.0 / 6075.0, 1804.0 / 7047.0, 1978.0 / 8091.0,
                    2350.0 / 9207.0, 2548.0 / 10395.0, 2968.0 / 11655.0, 3190.0 / 12987.0, 3658.0 / 14391.0};
    double sum1 = d[0];
    double delold = 1.0;
    double termold = d[0];
    int i = 1;
    double del;
    double term;
    while (1) {
        del = 1.0 / (1.0 + d[i + 1] * v * delold);
        term = termold * (del - 1.0);
        sum1 = sum1 + term;
        i++;
        delold = del;
        termold = term;
        if (i <= 20 || (fabs(termold) > 0.000001)) {
            break;
        }
    }
    return sum1;
}
