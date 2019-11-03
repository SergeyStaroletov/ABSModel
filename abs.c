//
//  abs.c
//  abs_demo
//
//  Created by Sergey Staroletov on 20/09/2019.
//  Copyright © 2019 Sergey Staroletov. All rights reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double t = 0;
double tMax = 20;
double dt = 0.01;

// constants
const double DesiredRelativeSlip = 0.2;
const double Pbmax = 1500;
const double Kf = 1;
const double I = 5;
const double v0 = 88;
const double Rr = 1.2500;
const double m = 50;
const double g = 32.180;
const double Ctrl = 1;
const double eps = 0.001;
const double TB = 0.01;

double y = 0;     // current value of the signal in the scheme
double yOld = 0;  // its old value
double btPart = 0;

// forward functions declarations
void tick(void);
double *mit(int m, double *x, double *y);
double lookUp(double x);
int sign(double x);
double HG(double x, double t);

int main(int argc, const char *argv[]) {
    //mi = mit(N, FX, FY); //for the spline implementation of lookup table
    
    t = 0;
    y = DesiredRelativeSlip;  // entry point
    
    while (t <= tMax) {
        tick();
        t = t + dt;
    }
    return 0;
}

// initial values for integrators
double int1_x = 0;
double int2_x = 0;
double int1_y = 0;
double int2_y = v0 / Rr;
double int3_x = 0;
double int3_y = v0;
double int4_x = 0;
double int4_y = 0;

void tick() {
    double ws;
    double vsAng;
    
    yOld = y;
    
    y = DesiredRelativeSlip;
    y -= yOld * Ctrl;
    
    // Bang-bang controller
    y = sign(y);
    
    // Hydraulic Lag
    y = HG(y, t);
    
    // Brake pressure
    double save = y;
    y = int1_y + (y + int1_x) * dt / 2;
    
    if (y < 0) y = 0;  // saturation
    if (y > Pbmax) y = Pbmax;
    int1_y = y;
    int1_x = save;
    
    // Force & torque
    y *= Kf;
    
    // Tire torque
    double part1 = lookUp(yOld);
    part1 *= m * g / 4;
    
    // Brake torque
    y = -y;
    y += part1 * Rr;
    
    y *= 1.0 / I;
    
    save = y;
    y = int2_y + (y + int2_x) * dt / 2;
    
    if (y < 0) y = 0;  // saturation
    if (y > 1000) y = 1000;
    int2_x = save;
    int2_y = y;
    
    ws = y;
    
    double y1 = y;
    double part2 = part1 * (-1) / m;
    
    save = part2;
    part2 = int3_y + (part2 + int3_x) * dt / 2;
    
    int3_x = save;
    
    if (part2 < 0) part2 = 0;  // saturation
    if (part2 > 1000) part2 = 1000;
    int3_y = part2;
    
    part2 *= 1 / Rr;
    
    vsAng = part2;
    
    double y2 = part2;
    double part3 = part2;
    part3 = part3 * dt;
    
    // double y3 = part3; //not used in the formulas
    
    // y = 1.0 - y1 / (y2 + (y2 == 0) * eps);
    y = 1.0 - y1 / (y2 + (fabs(y2) < eps) * eps);
    
    // printf("t: %lf   wheel speed: %lf    vehicle speed: %lf\n", t, ws, vsAng);
    printf("%lf;%lf;%lf\n", t, ws, vsAng);
}

int sign(double x) { return (x >= 0) ? 1 : -1; }

double f(double x, double y) {
    double k = 100;
    double T = TB;
    
    return (k * x - y) / T;
}

double rk_0 = 0;
double h = 0.001;

// runge kutta method
double RK(double x) {
    double y = rk_0;
    double k1 = f(x, y);
    double k2 = f(x + h / 2, y + h / 2 * k1);
    double k3 = f(x + h / 2, y + h / 2 * k2);
    double k4 = f(x + h, y + h * k3);
    
    double rk = rk_0 + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    rk_0 = rk;
    return rk;
}

double HG(double x, double t) {
    // 100 / (TBs + 1);
    // double k = 100;
    // double T = TB;
    // return x * k * (1 - exp(- t / T));
    return RK(x);
}

const double a = 0;
const double b = 1;
const double param = 0.3;

double FY[] = {0.0,   0.4,  0.8,  0.97, 1.0,  0.98, 0.96, 0.94, 0.92, 0.9, 0.88,
    0.855, 0.83, 0.81, 0.79, 0.77, 0.75, 0.73, 0.72, 0.71, 0.7};

double FX[] = {0.0,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5,
    0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95, 1.0};

// implementation of the lookup table as a constant step function
double lookupFake(double x) {
    double y = 0;
    if (x <= 0) y = 0;
    if (x <= 0.05) y = 0.4;
    if (x <= 0.1) y = 0.8;
    if (x <= 0.15) y = 0.97;
    if (x <= 0.2) y = 1.0;
    if (x <= 0.25) y = 0.98;
    if (x <= 0.3) y = 0.96;
    if (x <= 0.35) y = 0.94;
    if (x <= 0.4) y = 0.92;
    if (x <= 0.45) y = 0.90;
    if (x <= 0.50) y = 0.88;
    if (x <= 0.55) y = 0.855;
    if (x <= 0.60) y = 0.83;
    if (x <= 0.65) y = 0.81;
    if (x <= 0.70) y = 0.79;
    if (x <= 0.75) y = 0.77;
    if (x <= 0.80) y = 0.75;
    if (x <= 0.85) y = 0.73;
    if (x <= 0.90) y = 0.72;
    if (x <= 0.95) y = 0.71;
    if (x >= 1.0) y = 0.7;
    
    return y;
}

int N = sizeof(FY) / sizeof(double);
// Spline implementation for the lookup table (copy-pasted)
double *gaus(int n, double **a, double *b) {
    double *x, t, s;
    x = malloc(n * sizeof(double));
    int i, j, k;
    for (k = 0; k < n - 1; k++) {
        for (i = k + 1; i < n; i++) {
            t = a[i][k] / a[k][k];
            b[i] = b[i] - t * b[k];
            for (j = 0; j < n; j++) a[i][j] = a[i][j] - t * a[k][j];
        }
    }
    for (k = n - 1; k >= 0; k--) {
        s = 0;
        for (j = k + 1; j < n; j++) s += a[k][j] * x[j];
        x[k] = (b[k] - s) / a[k][k];
    }
    return x;
}

double proiz(int i, int m, double *x, double *y) {
    double res;
    if (i == 0)
        res = (y[1] - y[0]) / (x[1] - x[0]);
    else if (i == m - 1)
        res = (y[m - 1] - y[m - 2]) / (x[m - 1] - x[m - 2]);
    else
        res = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]);
    return res;
}

double *mit(int m, double *x, double *y) {
    int i, j;
    double *mi, **a, *b, h;
    a = malloc(m * sizeof(double *));
    for (i = 0; i < m; i++) a[i] = malloc(m * sizeof(double));
    b = malloc(m * sizeof(double));
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) a[i][j] = 0;
        if (i == 0 || i == m - 1)
            a[i][i] = 1;
        else {
            a[i][i - 1] = 1;
            a[i][i] = 4;
            a[i][i + 1] = 1;
            h = x[i] - x[i - 1];
            if (i != 0 && i != m - 1)
                b[i] =
                3 * (y[i + 1] - y[i - 1]) / h;
        }
        if (i == 0)
            b[i] = proiz(0, m, x, y);
        else if (i == m - 1)
            b[i] = proiz(m - 1, m, x, y);
    }
    mi = gaus(m, a, b);
    if (a) {
        for (i = 0; i < m; i++) free(a[i]);
        free(a);
    }
    if (b) free(b);
    return mi;
}

double spline_it(double xp, int m, double *x, double *y, double *mi) {
    int i;
    double sp = 0, drob1, drob2, drob3, drob4, xstep;
    
    for (i = 1; i < m; i++) {
        if (x[i - 1] <= xp && xp <= x[i]) {
            xstep = x[i] - x[i - 1];
            drob1 = (xp - x[i]) * (xp - x[i]) * (2 * (xp - x[i - 1]) + xstep) /
            pow(xstep, 3);
            drob2 = (xp - x[i - 1]) * (xp - x[i - 1]) * (2 * (x[i] - xp) + xstep) /
            pow(xstep, 3);
            drob3 = (xp - x[i]) * (xp - x[i]) * (xp - x[i - 1]) / pow(xstep, 2);
            drob4 = (xp - x[i - 1]) * (xp - x[i - 1]) * (xp - x[i]) / pow(xstep, 2);
            sp = y[i - 1] * drob1 + y[i] * drob2 + mi[i - 1] * drob3 + mi[i] * drob4;
        }
    }
    return sp;
}

double *mi;

double spline(double xp, int m, double *x, double *y) {
    double *mi, sp;
    mi = mit(m, x, y);
    sp = spline_it(xp, m, x, y, mi);
    // if(mi) delete [] mi;
    return sp;
}

double lookUp(double x) {
    if (x < FX[0]) return FY[0];
    if (x > FX[N - 1]) return FY[N - 1];
    // return spline_it(x, N, FX, FY, mi); //uncomment for spline implementation of
    // lookup table
    return lookupFake(x);  // we use simple solution to work as a correspandant
    // solution in KeY syntanx
    
    /* //Lagrange implementation of interplataion is also commented
     double rez = 0;
     for (int k = 0; k < N; k++) {
     double rezz = 1;
     for (int n  = 0; n < N; n++)
     if (n != k)
     rezz *= (x - FX[n]) / (FX[k] - FX[n]);
     rez += FY[k] * rezz;
     }
     rez = fabs(rez);
     printf("look(%lf) = %lf\n", x, rez);
     return rez;
     */
}
