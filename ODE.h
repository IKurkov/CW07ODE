/* Kurkov Ivan, 22.B05-MM, 12.05.2024 */
#ifndef ODE_H
#define ODE_H

struct ODEUnit
{
  const char *Name;
  double (*Method)(double, double, double, double(*)(double, double));
};

double Euler( double x_k, double y_k, double h, double (*f)(double, double) );
double EulerMidRecMod( double x_k, double y_k, double h, double (*f)(double, double) );
double EulerTrapezeMod( double x_k, double y_k, double h, double (*f)(double, double) );
double RungeKutta4( double x_n, double y_n, double h, double (*f)(double, double) );
double ExtAdams4( double x_end, const double *y_end, double h, double (*f)(double, double) );

extern const ODEUnit ODEList[];
extern const size_t ODElen;

#endif // !ODE_H

