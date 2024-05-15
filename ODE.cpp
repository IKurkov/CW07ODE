/* Kurkov Ivan, 22.B05-MM, 12.05.2024 */
#include "ode.h"

const ODEUnit ODEList[] =
{
  {"Euler", Euler},
  {"Euler [ * ]", EulerMidRecMod},
  {"Euler trapeze", EulerTrapezeMod},
  {"Runge-Kutta 4", RungeKutta4}
};
const size_t ODElen = 4;

/* Euler's method for solving Cauchy problem:
   y_{k + 1} = y_k + h * f(x_k, y_k)*/
double Euler( double x_k, double y_k, double h, double (*f)(double, double) )
{
  return y_k + h * f(x_k, y_k);
}

/* Euler's method middle rectangles modification for solving Cauchy problem:
   y_{k + 1/2} = y_k + h / 2 * f(x_k, y_k)
   y_{k + 1} = y_k + h * f(x_k + h / 2, y_{k + 1/2})*/
double EulerMidRecMod( double x_k, double y_k, double h, double (*f)(double, double) )
{
  return y_k + h * f(x_k + h / 2, y_k + h * f(x_k, y_k) / 2);
}

/* Euler's method trapeze modification for solving Cauchy problem:
   u_{k + 1} = y_k + h * f(x_k, y_k)
   y_{k + 1} = y_k + h / 2 * (f(x_k, y_k) + f(x_{k + 1}, u_{k + 1}))*/
double EulerTrapezeMod( double x_k, double y_k, double h, double (*f)(double, double) )
{
  double field = f(x_k, y_k);

  return y_k + h / 2 * (f(x_k, y_k) + f(x_k + h, y_k + h * field));
}

/* Runge-Kutta method for solving Cauchy problem */
double RungeKutta4( double x_n, double y_n, double h, double (*f)(double, double) )
{
  double k1 = h * f(x_n, y_n),
    k2 = h * f(x_n + h / 2, y_n + k1 / 2),
    k3 = h * f(x_n + h / 2, y_n + k2 / 2),
    k4 = h * f(x_n + h, y_n + k3);

  return y_n + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

/* Extrapolation Adams method for solving Cauchy problem */
double ExtAdams4( double x_end, const double *y_end, double h, double (*f)(double, double) )
{
  double finite_dif[5];
  
  for (int i = 4; i >= 0; i--, x_end -= h)
    finite_dif[i] = h * f(x_end, y_end[i]);
  for (int i = 4; i >= 1; i--)
    for (int j = 0; j < i; j++)
      finite_dif[j] = finite_dif[j + 1] - finite_dif[j];
  return y_end[4] + finite_dif[4] + finite_dif[3] / 2 + 5 * finite_dif[2] / 12
    + 3 * finite_dif[1] / 8 + 251 * finite_dif[0] / 720;
}