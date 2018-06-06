#include <cmath>

#include "global_definitions.h"
#include "Fluid1D.h"

Real Sign(Real x) {
  if(x < 0.0){
    return -1.0;
  } else {
    return 1.0;
  }
}

Real Min(Real a, Real b, Real c){
  return std::fmin(std::fmin(a, b), c);
}

Real Gradient(const Rvec &x, const Rvec &sia, const Rvec &via, Int32 i)
{
  Real dx[3];
  Real f[3];

  for(auto j = 0; j < 3; ++j){
    auto iv = i - 1+ j;
    dx[j] = x[iv + 1] - x[iv];
    f[j] = 1.5 * via[iv] - 0.25 * (sia[iv + 1] + sia[iv]);
  }

  auto f21 = f[2] - f[1], f10 = f[1] - f[0];
  auto df = (dx[1] / (dx[0] + dx[1] + dx[2]))
          * ((2.0 * dx[0] + dx[1]) / (dx[1] + dx[2]) * f21
          + (dx[1] + 2.0 * dx[2]) / (dx[0] + dx[1]) * f10);
  auto cond = f21 * f10 <= 0.0;
  if (cond) {
    return 0.0;
  } else {
    return Sign(df) / dx[1]
                    * Min(std::fabs(df),
                          2.0 * std::fabs(f21),
                          2.0 * std::fabs(f10));
  }
}

Real Fluid1D::CipCsl3(Int32 i, Real v,
                      const Rvec &sia, const Rvec &via) noexcept
{
  auto xi = -v * dt_;

  Real c[4];
  Int32 ip, iv;
  Real dx;
  if(v < 0.0){
    ip = i + 1;
    iv = i;
    dx = x_[ip] - x_[i];
  }else{
    ip = i - 1;
    iv = i - 1;
    dx = x_[ip] - x_[i];
  }
  auto d = Gradient(x_, sia, via, iv);
  c[0] = sia[i];
  c[1] = 2.0 / dx * (3.0 * via[iv] - 3.0 * sia[i] - d * dx);
  c[2] = 3.0 / (dx * dx)
       * (-2.0 * via[iv] + 3.0 * sia[i] - sia[ip] + 2.0 * d * dx);
  c[3] = 4.0 / (dx * dx * dx) * (sia[ip] - sia[i] - d * dx);

  return c[0] + c[1] * xi + c[2] * xi * xi + c[3] * xi * xi * xi;
}

Real Fluid1D::CipCsl3(Int32 i, Real v, Real v_bck, Real v_fwd,
                      const Rvec &sia, const Rvec &via) noexcept
{
  auto dx_fwd = x_[i + 1] - x_[i];
  auto dx_bck = x_[i] - x_[i - 1];
  auto du_dx = (v_fwd - v_bck) / (dx_fwd + dx_bck);
  auto d2u_dx2 = 2.0 * (v_fwd - 2.0 * v + v_bck)
               / (dx_fwd * dx_fwd + dx_bck * dx_bck);

  auto xi = -v * dt_ + 0.5 * v * (dt_ * dt_) * du_dx
          - v * std::pow(dt_, 3) / 6.0
          * (du_dx * du_dx + v * d2u_dx2);

  Real c[4];
  Int32 ip, iv;
  Real dx;
  if(xi > 0.0){
    ip = i + 1;
    iv = i;
    dx = x_[ip] - x_[i];
  }else{
    ip = i - 1;
    iv = i - 1;
    dx = x_[ip] - x_[i];
  }
  auto d = Gradient(x_, sia, via, iv);
  c[0] = sia[i];
  c[1] = 2.0 / dx * (3.0 * via[iv] - 3.0 * sia[i] - d * dx);
  c[2] = 3.0 / (dx * dx)
       * (-2.0 * via[iv] + 3.0 * sia[i] - sia[ip] + 2.0 * d * dx);
  c[3] = 4.0 / (dx * dx * dx) * (sia[ip] - sia[i] - d * dx);

  return c[0] + c[1] * xi + c[2] * xi * xi + c[3] * xi * xi * xi;
}
