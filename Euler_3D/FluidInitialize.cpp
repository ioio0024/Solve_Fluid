#include <cmath>

#include "global_definitions.h"
#include "Vector3D.h"
#include "FluidInitialize.h"

void FluidInitialize::operator ()(const Rvec &x, const Rvec &y, const Rvec &z,
                                  Vector3DR *d,
                                  Vector3DR *u,
                                  Vector3DR *v,
                                  Vector3DR *w,
                                  Vector3DR *p) noexcept
{
  Int32 nx = x.size();
  Int32 ny = y.size();
  Int32 nz = z.size();
  auto x_avg = 0.5 * (x[GHOST] + x[nx - GHOST - 1]);
  auto y_avg = 0.5 * (y[GHOST] + y[ny - GHOST - 1]);
  auto z_avg = 0.5 * (z[GHOST] + z[nz - GHOST - 1]);
  auto r_max = std::sqrt(std::pow(x[GHOST] - x_avg, 2)
                         + std::pow(y[GHOST] - y_avg, 2)
                         + std::pow(z[GHOST] - z_avg, 2));

  for(auto k = 0; k < nz; ++k){
    for(auto j = 0; j < ny; ++j){
      for(auto i = 0; i < nx; ++i){
        auto r = std::sqrt(std::pow(x[i] - x_avg, 2)
                         + std::pow(y[j] - y_avg, 2)
                         + std::pow(z[k] - z_avg, 2));
        if(r <= 0.2 * r_max) {
          (*d)(i, j, k) = 1.0;
          (*p)(i, j, k) = 1.0;
        } else {
          (*d)(i, j, k) = 0.125;
          (*p)(i, j, k) = 0.1;
        }
        (*u)(i, j, k) = 0.0;
        (*v)(i, j, k) = 0.0;
        (*w)(i, j, k) = 0.0;
      }
    }
  }
}
