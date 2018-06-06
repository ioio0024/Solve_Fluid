#include <cmath>

#include "global_definitions.h"
#include "Vector2D.h"
#include "FluidInitialize.h"

void FluidInitialize::operator ()(const Rvec &x, const Rvec &y,
                                  Vector2D *d,
                                  Vector2D *u,
                                  Vector2D *v,
                                  Vector2D *p) noexcept
{
  Int32 nx = x.size();
  Int32 ny = y.size();
  auto x_avg = 0.5 * (x[GHOST] + x[nx - GHOST - 1]);
  auto y_avg = 0.5 * (y[GHOST] + y[ny - GHOST - 1]);
  auto r_max = std::sqrt(std::pow(x[GHOST] - x_avg, 2)
                         + std::pow(y[GHOST] - y_avg, 2));

  for(auto j = 0; j < ny; ++j){
    for(auto i = 0; i < nx; ++i){
      auto r = std::sqrt(std::pow(x[i] - x_avg, 2) + std::pow(y[j] - y_avg, 2));
      if(r <= 0.2 * r_max) {
        (*d)(i, j) = 1.0;
        (*p)(i, j) = 1.0;
      } else {
        (*d)(i, j) = 0.125;
        (*p)(i, j) = 0.1;
      }
      (*u)(i, j) = 0.0;
      (*v)(i, j) = 0.0;
    }
  }
}
