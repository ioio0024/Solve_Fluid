#include <cmath>

#include "FluidInitialize.h"

void FluidInitialize::operator ()(
    const Vector3DR &x, const Vector3DR &y, const Vector3DR &z,
    FluidQuantity3DArray *q) noexcept
{
  Int32 nx = x.size();
  Int32 ny = y.size();
  Int32 nz = z.size();
  auto x_avg = 0.5 * (x(GHOST) + x(nx - GHOST - 1));
  auto y_avg = 0.5 * (y(GHOST) + y(ny - GHOST - 1));
  auto z_avg = 0.5 * (z(GHOST) + z(nz - GHOST - 1));
  auto r_max = std::sqrt(std::pow(x(GHOST) - x_avg, 2)
                         + std::pow(y(GHOST) - y_avg, 2)
                         + std::pow(z(GHOST) - z_avg, 2));

  for(auto k = 0; k < nz; ++k){
    for(auto j = 0; j < ny; ++j){
      for(auto i = 0; i < nx; ++i){
        auto r = std::sqrt(std::pow(x(i) - x_avg, 2)
                         + std::pow(y(j) - y_avg, 2)
                         + std::pow(z(k) - z_avg, 2));
        if(r <= 0.2 * r_max) {
          q->d(i, j, k) = 1.0;
          q->p(i, j, k) = 1.0;
        } else {
          q->d(i, j, k) = 0.125;
          q->p(i, j, k) = 0.1;
        }
        q->u(i, j, k) = 0.0;
        q->v(i, j, k) = 0.0;
        q->w(i, j, k) = 0.0;
      }
    }
  }
}
