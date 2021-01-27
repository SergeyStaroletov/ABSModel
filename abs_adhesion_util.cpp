/*
 * ABS model testing demo according to EU regulations
 * (c) Sergey Staroletov, 2021
 */
#include <cmath>
#include <iostream>
using namespace std;

/*
 * Run ABS braking  cycle: from initial speed v0
 * with or without abs operation
 * control time from vControl1 to vControl2
 * return time of stop and control time and the distance
 */
void ABS(double v0, bool abs, double vControl1, double vControl2,
         double &retTStop, double &retTControl, double &retSStop);

// For RK
double rk_0 = 0;

// Vehicle constants
const double m = 1355 + 75;
const double g = 9.8;
double P = m;
double H = 0.5;
double E = 2.675;
double Fr = 0.4 * m;
double Ff = 0.6 * m;

int main(__unused int argc, __unused const char *argv[]) {
  double retTStop = 0;
  double retTControl = 0;
  double retSStop = 0;

  /*
     Given an initial vehicle speed of 55 km/h, the maximum braking rate
     ($z_{AL}$) is measured with full cycle of the anti-lock
     braking system using the time taken $t_{m_1}$  to reduce the speed from 45
     km/h to 15 km/h
  */

  double vControl1 = 45;
  double vControl2 = 15;
  double v0 = 55;
  // Run our ABS model
  ABS(v0, true, vControl1, vControl2, retTStop, retTControl, retSStop);

  double tm1 = retTControl;
  cout << "tm1 = " << tm1 << endl;

  double zal = 0.849 / tm1;
  cout << "z_al = " << zal << endl;

  /*
     Dynamic axle loads are calculated on static loads $F_f$ and $F_r$ with
     respect to gravity center $h$ and wheel base $E$:
  */
  double Frdyn = Fr - H / E * zal * P * g;
  double Ffdyn = Ff + H / E * zal * P * g;

  cout << "Frdyn = " << Frdyn << endl;
  cout << "Ffdyn = " << Ffdyn << endl;

  /*
     The brakes without locking the wheels shall be applied on only one axle of
     the vehicle under special tests, at an initial speed of 50 km/h A number of
     tests at increments of line pressure shall be carried out to determine the
     maximum braking rate of the vehicle ($z_m$). During each test, a constant
     input force shall be maintained to determine the braking rate by reference
     to the minimal time taken ($t_{m_2}$) for the speed to reduce from 40 km/h
     to 20 km/h using the formula:  z_m = \dfrac {0,566}{t_{m_2}}
  */
  vControl1 = 40;
  vControl2 = 20;
  v0 = 50;
  // The anti-lock system in the following tests shall be turned off
  ABS(v0, false, vControl1, vControl2, retTStop, retTControl, retSStop);

  double tm2 = retTControl;
  cout << "tm2 = " << tm2 << endl;

  double zm = 0.566 / tm2;
  cout << "zm = " << zm << endl;

  double kf = (zm * P * g - 0.015 * Fr) / (Ff + H / E * zm * P * g);
  double kr = (zm * P * g - 0.010 * Ff) / (Fr - H / E * zm * P * g);

  cout << "kf = " << kf << endl;
  cout << "kr = " << kr << endl;

  /*
     The coefficient of adhesion $k_M$  is determined by taking into
     account the dynamic axle loads (front and rear): k_M = \dfrac{k_f \cdot
     F_{f_{dyn}}
     + k_r \cdot F_{r_{dyn}}}{P \cdot g}
  */
  double kM = (kf * Ffdyn + kr * Frdyn) / (P * g);

  cout << "kM = " << kM << endl;

  /*
     With respect to the document, $\epsilon$ is defined as the quotient of the
     maximum braking rate with the anti-lock system operative ($z_{AL}$) and the
     coefficient of adhesion ($k_M$): \epsilon = \dfrac{z_{AL}}{k_M}
  */

  double epsilon = zal / kM;
  cout << "Epsilon = " << epsilon << endl;

  if (epsilon < 0.75)
    cout << "ABS works ineffectively" << endl;
  else
    cout << "ABS works effectively" << endl;

  return 0;
}

// Forward functions declarations
double lookUp(double x);
int sign(double x);
double HG(double x, double t);

void ABS(double v0, bool abs, double vControl1, double vControl2,
         double &retTStop, double &retTControl, double &retSStop) {
  double t = 0;
  double tMax = 20;
  double dt = 0.001;
  double msToKmh = 1000.0 / 3600.0;
  v0 *= msToKmh;
  vControl1 *= msToKmh;
  vControl2 *= msToKmh;

  // Constants for the model
  const double DesiredRelativeSlip = 0.2;
  const double Pbmax = 1500;
  const double Kf = 2.5;
  const double I = 5;
  const double Rr = 0.35;
  const double Ctrl = abs ? 1 : 0;
  const double eps = 0.001;
  const double Kwet = 1;

  double y = 0;     // current value of the signal in the scheme
  double yOld = 0;  // its old value

  // Initial values for integrators
  double int1_x = 0;
  double int2_x = 0;
  double int1_y = 0;
  double int2_y = v0 / Rr;
  double int3_x = 0;
  double int3_y = v0;
  rk_0 = 0;

  double t1 = 0;
  double t2 = 0;

  double ws = v0 / Rr;
  double vsAng;
  double v = v0;
  t = 0;
  y = DesiredRelativeSlip;  // entry point

  // Time loop
  while (t <= tMax && v >= eps) {
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
    double part1 = lookUp(yOld) / Kwet;
    part1 *= m * g / 4;

    // Brake torque
    y = -y;
    y += part1 * Rr;
    y *= 1.0 / I;

    save = y;
    y = int2_y + (y + int2_x) * dt / 2;

    if (y < 0) y = 0;
    if (y > 1000) y = 1000;
    int2_x = save;
    int2_y = y;

    ws = y;

    double y1 = y;
    double part2 = part1 * (-1) / m;

    save = part2;
    part2 = int3_y + (part2 + int3_x) * dt / 2;
    int3_x = save;

    if (part2 < 0) part2 = 0;
    if (part2 > 1000) part2 = 1000;
    int3_y = part2;

    part2 *= 1 / Rr;
    vsAng = part2;

    double y2 = part2;
    double part3 = part2;
    part3 = part3 * dt;

    retSStop = part3;

    y = 1.0 - y1 / (y2 + (fabs(y2) < eps) * eps);

    v = vsAng * Rr;

    // Check conditions to save the time
    if (fabs(v - vControl1) < 0.5) t1 = t;
    if (fabs(v - vControl2) < 0.5) t2 = t;

    cout << t << "," << ws * Rr / msToKmh << "," << v / msToKmh << endl;

    t = t + dt;
  }
  retTStop = t;
  retTControl = (t2 - t1) / 2;
}

int sign(double x) { return (x >= 0) ? 1 : -1; }

// For RK
double f(double x, double y) {
  double k = 100;
  double T = 0.05;

  return (k * x - y) / T;
}

const double h = 0.001;

// Runge Kutta method
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

double HG(double x, __unused double t) {
  // 100 / (TBs + 1);
  // double k = 100;
  // double T = TB;
  // return x * k * (1 - exp(- t / T));
  return RK(x);
}

double FY[] = {0.0,   0.4,  0.8,  0.97, 1.0,  0.98, 0.96, 0.94, 0.92, 0.9, 0.88,
               0.855, 0.83, 0.81, 0.79, 0.77, 0.75, 0.73, 0.72, 0.71, 0.7};

double FX[] = {0.0,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5,
               0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95, 1.0};

// Implementation of the lookup table as a constant step function
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

const int N = sizeof(FY) / sizeof(double);

double lookUp(double x) {
  if (x < FX[0]) return FY[0];
  if (x > FX[N - 1]) return FY[N - 1];
  // return spline_it(x, N, FX, FY, mi); //uncomment for spline implementation
  // or lookup table
  return lookupFake(x);
}
