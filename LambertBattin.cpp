
//
// Created by josej on 24/03/2023.
//

#include "LambertBattin.h"
#include <cmath>
#include <cstring>
#include "seebattk.h"
#include "seebatt.h"

void LambertBattin(const std::vector<double>& ro, const std::vector<double>& r, const char* dm, double Dtsec,
                   std::vector<double>& vo, std::vector<double>& v) {
    const double small = 0.000001;
    const double mu = 3.986004418e14;
    const double magr = norm(r, 3);
    const double magro = norm(ro, 3);
    const double CosDeltaNu = dot(ro, r, 3, 3) / (magro * magr);
    const std::vector<double> rcrossr = cross(ro, r, 3, 3);
    const double magrcrossr = norm(rcrossr, 3);
    double SinDeltaNu;
    double DNu;
    double Ror;
    double eps;
    double tan2w;
    double rp;
    double L;
    double m;
    double x;
    double xn;
    double chord;
    double s;
    double lim1;
    int Loops;
    double tempx;
    double Denom;
    double h1;
    double h2;
    double b;
    double u;
    double k2;
    double y;
    double a;
    double arg1;
    double arg2;
    double AlpH;
    double BetH;
    double DH;
    double F;
    double GDot;
    double G;
    double Sinv;
    double Cosv;
    double BetE;
    double am;
    double ae;
    double be;
    double tm;
    double AlpE;
    double DE;
    if (std::strcmp(dm, "pro") == 0) {
        SinDeltaNu = magrcrossr / (magro * magr);
    } else {
        SinDeltaNu = -magrcrossr / (magro * magr);
    }
    DNu = std::atan2(SinDeltaNu, CosDeltaNu);
    if (DNu < 0.0) {
        DNu += 2.0 * M_PI;
    }
    Ror = magr / magro;
    eps = Ror - 1.0;
    tan2w = 0.25 * eps * eps / (std::sqrt(Ror) + Ror * (2.0 + std::sqrt(Ror)));
    rp = std::sqrt(magro * magr) * (std::pow(std::cos(DNu * 0.25), 2) + tan2w);
    if (DNu < M_PI) {
        L = (std::pow(std::sin(DNu * 0.25), 2) + tan2w) / ((std::pow(std::sin(DNu * 0.25), 2)) + tan2w + std::cos(DNu * 0.5));
    } else {
        L = (std::pow(std::cos(DNu * 0.25), 2) + tan2w - std::cos(DNu * 0.5)) / ((std::pow(std::cos(DNu * 0.25), 2)) + tan2w);
    }
    m = mu * Dtsec * Dtsec / (8.0 * rp * rp * rp);
    x = 10.0;
    xn = L;
    chord = std::sqrt(magro * magro + magr * magr - 2.0 * magro * magr * std::cos(DNu));
    s = (magro + magr + chord) * 0.5;
    lim1 = std::sqrt(m / L);
    Loops = 1;
    while (true) {
        x = xn;
        tempx = seebatt(x);
        Denom = 1.0 / ((1.0 + 2.0 * x + L) * (4.0 * x + tempx * (3.0 + x)));
        h1 = std::pow(L + x, 2) * (1.0 + 3.0 * x + tempx) * Denom;
        h2 = m * (x - L + tempx) * Denom;
        b = 0.25 * 27.0 * h2 / (std::pow(1.0 + h1, 3));
        if (b < -1.0) {
            xn = 1.0 - 2.0 * L;
        } else {
            if (y1 > lim1) {
                xn *= (lim1 / y1);
            } else {
                u = 0.5 * b / (1.0 + std::sqrt(1.0 + b));
                k2 = seebattk(u);
                y = ((1.0 + h1) / 3.0) * (2.0 + std::sqrt(1.0 + b) / (1.0 + 2.0 * u * k2 * k2));
                xn = std::sqrt(std::pow((1.0 - L) * 0.5, 2) + m / (y * y)) - (1.0 + L) * 0.5;
            }
        }
        Loops++;
        y1 = std::sqrt(m / ((L + x) * (1.0 + x)));
        if ((std::fabs(xn - x) < small) && (Loops > 30)) {
            break;
        }
    }
    a = mu * Dtsec * Dtsec / (16.0 * rp * rp * xn * y * y);
    if (a < -small) {
        arg1 = std::sqrt(s / (-2.0 * a));
        arg2 = std::sqrt((s - chord) / (-2.0 * a));
        AlpH = 2.0 * std::asinh(arg1);
        BetH = 2.0 * std::asinh(arg2);
        DH = AlpH - BetH;
        F = 1.0 - (a / magro) * (1.0 - std::cosh(DH));
        GDot = 1.0 - (a / magr) * (1.0 - std::cosh(DH));
        G = Dtsec - std::sqrt(-a * a * a / mu) * (std::sinh(DH) - DH);
    } else {
        if (a > small) {
            arg1 = std::sqrt(s / (2.0 * a));
            arg2 = std::sqrt((s - chord) / (2.0 * a));
            Sinv = arg2;
            Cosv = std::sqrt(1.0 - (magro + magr - chord) / (4.0 * a));
            BetE = 2.0 * std::asin(Cosv);
            BetE = 2.0 * std::asin(Sinv);
            if (DNu > M_PI) {
                BetE = -BetE;
            }
            Cosv = std::sqrt(1.0 - s / (2.0 * a));
            Sinv = arg1;
            am = s * 0.5;
            ae = M_PI;
            be = 2.0 * std::asin(std::sqrt((s - chord) / s));
            tm = std::sqrt(am * am * am / mu) * (ae - (be - std::sin(be)));
            if (Dtsec > tm) {
                AlpE = 2.0 * M_PI - 2.0 * std::asin(Sinv);
            } else {
                AlpE = 2.0 * std::asin(Sinv);
            }
            DE = AlpE - BetE;
            F = 1.0 - (a / magro) * (1.0 - std::cos(DE));
            GDot = 1.0 - (a / magr) * (1.0 - std::cos(DE));
            G = Dtsec - std::sqrt(a * a * a / mu) * (DE - std::sin(DE));
        } else {
            arg1 = 0.0;
            arg2 = 0.0;
            throw std::runtime_error("a parbolic orbit");
        }
    }
    for (int i = 0; i < 3; i++) {
        vo[i] = (r[i] - F * ro[i]) / G;
        v[i] = (GDot * r[i] - ro[i]) / G;
    }
}
