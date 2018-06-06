#include "FluidInitialize.h"

void FluidInitialize(const Rvec &x,
                     Rvec *density, Rvec *velocity, Rvec *pressure)
{
  auto n = x.size();
  auto xmin = x[GHOST], xmax = x[n - GHOST - 1];
  auto x_avg = (xmax + xmin) / 2.0;
  for(auto i = 0; i < n; ++i){
    Real rho, u, p;
    if(x[i] <= x_avg){
      rho = 1.0;
      p = 1.0;
    }else{
      rho = 0.125;
      p = 0.1;
    }
    u = 0.0;

    (*density)[i] = rho;
    (*velocity)[i] = u;
    (*pressure)[i] = p;
  }
}
