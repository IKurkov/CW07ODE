/* Kurkov Ivan, 22.B05-MM, 12.05.2024 */
#include <cmath>
#include <conio.h>
#include <iostream>

#include "fort.hpp"
#include "ode.h"

const char *V7Str = "y' = -y + cos(x)";

double V7ODE( double x, double y )
{
  return -y + cos(x);
}

double V7Analytic( double x, double x0, double y0 )
{
  return (y0 - 0.5 * (sin(x0) + cos(x0))) * exp(x0 - x) + 0.5 * (sin(x) + cos(x));
}

double V7Taylor( double x, double x0, double y0 )
{
  return y0 + (cos(x0) - y0) * x + (y0 - cos(x0) - sin(x0)) * pow(x, 2) / 2 + (sin(x0) - y0) * pow(x, 3) / 6
    + y0 * pow(x, 4) / 24 + (cos(x0) - y0) * pow(x, 5) / 120 + (y0 - cos(x0) - sin(x0)) * pow(x, 6) / 720;
}

int main( void )
{
  size_t N;
  double h, x0 = 0, y0 = 1, x, y, y_analytic, *y_vals = nullptr;
  bool run = true;

  while (run)
  {
    std::cout << "Numeric solution for Cauchy problem menu [variant #7]:\n"
      "0 - exit\n"
      "1 - solve Cauchy problem\n"
      "2 - set initial data\n";
    switch (_getch())
    {
    case '0':
      run = false;
      break;
    case '1':
    {
      fort::char_table Taylor, comparison;

      do
      {
        std::cout << "Input N >= 3: ";
        std::cin >> N;
      } while (N < 3);
      do
      {
        std::cout << "Input h > 0: ";
        std::cin >> h;
      } while (h <= 0);

      Taylor << fort::header << "x_i" << "y(x_i)" << "y_i[Taylor]" << "|y(x_N) - y_N|" << fort::endr;
      y_vals = new double[N + 3];
      x = x0 - 2 * h;
      for (size_t i = 0; i < 5; i++, x += h)
      {
        y_analytic = V7Analytic(x, x0, y0);
        y_vals[i] = y = V7Taylor(x, x0, y0);
        Taylor << x << std::setprecision(15) << y_analytic << y << fabs(y_analytic - y) << fort::endr;
      }
      for (size_t i = 5; i <= N + 2; i++, x += h)
      {
        y_analytic = V7Analytic(x, x0, y0);
        y = V7Taylor(x, x0, y0);
        Taylor << x << y_analytic << y << fabs(y_analytic - y) << fort::endr;
      }
      std::cout << "\nODE:" << V7Str << '\n' << Taylor.to_string();

      comparison << fort::header << "Method" << "y_N[Method]" << "|y(x_N) - y_N|" << fort::endr;
      comparison << "Taylor" << std::setprecision(15) << y << fabs(y_analytic - y) << fort::endr;
      for (size_t i = 0; i < ODElen; i++)
      {
        y = y0;
        x = x0;
        for (size_t j = 1; j <= N; j++, x += h)
          y = ODEList[i].Method(x, y, h, V7ODE);
        comparison << ODEList[i].Name << std::setprecision(15) << y << fabs(y_analytic - y) << fort::endr;
      }
      x = x0 + 2 * h;
      for (size_t i = 5; i <= N + 2; i++, x += h)
        y_vals[i] = ExtAdams4(x, y_vals + i - 5, h, V7ODE);
      comparison << "Extrapolation Adams 4" << std::setprecision(15) << y_vals[N + 1] << fabs(y_analytic - y_vals[N + 2]) << fort::endr;
      std::cout << comparison.to_string();

      delete[] y_vals;
      y_vals = nullptr;
      break;
    }
    case '2':
      std::cout << "ODE: " << V7Str << '\n';
      std::cout << "Input x0 (current value = " << x0 << "): ";
      std::cin >> x0;
      std::cout << "Input y0 (current value = " << y0 << "): ";
      std::cin >> y0;
      break;
    default:
      std::cout << "[Error]: Incorrect choice!\n";
      break;
    }
  }
  return 0;
}