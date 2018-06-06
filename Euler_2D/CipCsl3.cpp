#include <cassert>
#include <cmath>

#include "global_definitions.h"
#include "Fluid2D.h"

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

Real GradientX(const Rvec &x, const Vector2D &sia, const Vector2D &via,
               Int32 i, Int32 j)
{
  Real dx[3];
  Real f[3];

  for(auto n = 0; n < 3; ++n){
    auto iv = i - 1+ n;
    dx[n] = x[iv + 1] - x[iv];
    f[n] = 1.5 * via(iv, j) - 0.25 * (sia(iv + 1, j) + sia(iv, j));
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

Real GradientY(const Rvec &x, const Vector2D &sia, const Vector2D &via,
               Int32 i, Int32 j)
{
  Real dx[3];
  Real f[3];

  for(auto n = 0; n < 3; ++n){
    auto jv = j - 1+ n;
    dx[n] = x[jv + 1] - x[jv];
    f[n] = 1.5 * via(i, jv) - 0.25 * (sia(i, jv + 1) + sia(i, jv));
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

Real Fluid2D::CipCsl3X(Int32 i, Int32 j, Real v,
                       const Vector2D &sia, const Vector2D &via) noexcept
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
  assert(std::fabs(xi) <= std::fabs(dx));
  auto d = GradientX(x_, sia, via, iv, j);
  c[0] = sia(i, j);
  c[1] = 2.0 / dx * (3.0 * via(iv, j) - 3.0 * sia(i, j) - d * dx);
  c[2] = 3.0 / (dx * dx)
       * (-2.0 * via(iv, j) + 3.0 * sia(i, j) - sia(ip, j) + 2.0 * d * dx);
  c[3] = 4.0 / (dx * dx * dx) * (sia(ip, j) - sia(i, j) - d * dx);

  return c[0] + c[1] * xi + c[2] * xi * xi + c[3] * xi * xi * xi;
}

Real Fluid2D::CipCsl3Y(Int32 i, Int32 j, Real v,
                       const Vector2D &sia, const Vector2D &via) noexcept
{
  auto xi = -v * dt_;

  Real c[4];
  Int32 jp, jv;
  Real dy;
  if(v < 0.0){
    jp = j + 1;
    jv = j;
    dy = y_[jp] - y_[j];
  }else{
    jp = j - 1;
    jv = j - 1;
    dy = y_[jp] - y_[j];
  }
  assert(std::fabs(xi) <= std::fabs(dy));
  auto d = GradientY(y_, sia, via, i, jv);
  c[0] = sia(i, j);
  c[1] = 2.0 / dy * (3.0 * via(i, jv) - 3.0 * sia(i, j) - d * dy);
  c[2] = 3.0 / (dy * dy)
       * (-2.0 * via(i, jv) + 3.0 * sia(i, j) - sia(i, jp) + 2.0 * d * dy);
  c[3] = 4.0 / (dy * dy * dy) * (sia(i, jp) - sia(i, j) - d * dy);

  return c[0] + c[1] * xi + c[2] * xi * xi + c[3] * xi * xi * xi;
}
