#include <cassert>
#include <cmath>

#include "global_definitions.h"
#include "Fluid3D.h"

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

Real GradientX(const Rvec &x, const Vector3DR &sia, const Vector3DR &via,
               Int32 i, Int32 j, Int32 k)
{
  Real dx[3];
  Real f[3];

  for(auto n = 0; n < 3; ++n){
    auto iv = i - 1+ n;
    dx[n] = x[iv + 1] - x[iv];
    f[n] = 1.5 * via(iv, j, k) - 0.25 * (sia(iv + 1, j, k) + sia(iv, j, k));
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

Real GradientY(const Rvec &x, const Vector3DR &sia, const Vector3DR &via,
               Int32 i, Int32 j, Int32 k)
{
  Real dx[3];
  Real f[3];

  for(auto n = 0; n < 3; ++n){
    auto jv = j - 1+ n;
    dx[n] = x[jv + 1] - x[jv];
    f[n] = 1.5 * via(i, jv, k) - 0.25 * (sia(i, jv + 1, k) + sia(i, jv, k));
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

Real GradientZ(const Rvec &x, const Vector3DR &sia, const Vector3DR &via,
               Int32 i, Int32 j, Int32 k)
{
  Real dx[3];
  Real f[3];

  for(auto n = 0; n < 3; ++n){
    auto kv = k - 1+ n;
    dx[n] = x[kv + 1] - x[kv];
    f[n] = 1.5 * via(i, j, kv) - 0.25 * (sia(i, j, kv + 1) + sia(i, j, kv));
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

Real Fluid3D::CipCsl3X(Int32 i, Int32 j, Int32 k, Real v,
                       const Vector3DR &sia, const Vector3DR &via) noexcept
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
  auto d = GradientX(x_, sia, via, iv, j, k);
  c[0] = sia(i, j, k);
  c[1] = 2.0 / dx * (3.0 * via(iv, j, k) - 3.0 * sia(i, j, k) - d * dx);
  c[2] = 3.0 / (dx * dx)
       * (-2.0 * via(iv, j, k) + 3.0 * sia(i, j, k) - sia(ip, j, k) + 2.0 * d * dx);
  c[3] = 4.0 / (dx * dx * dx) * (sia(ip, j, k) - sia(i, j, k) - d * dx);

  return c[0] + c[1] * xi + c[2] * xi * xi + c[3] * xi * xi * xi;
}

Real Fluid3D::CipCsl3Y(Int32 i, Int32 j, Int32 k, Real v,
                       const Vector3DR &sia, const Vector3DR &via) noexcept
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
  auto d = GradientY(y_, sia, via, i, jv, k);
  c[0] = sia(i, j, k);
  c[1] = 2.0 / dy * (3.0 * via(i, jv, k) - 3.0 * sia(i, j, k) - d * dy);
  c[2] = 3.0 / (dy * dy)
       * (-2.0 * via(i, jv, k) + 3.0 * sia(i, j, k) - sia(i, jp, k) + 2.0 * d * dy);
  c[3] = 4.0 / (dy * dy * dy) * (sia(i, jp, k) - sia(i, j, k) - d * dy);

  return c[0] + c[1] * xi + c[2] * xi * xi + c[3] * xi * xi * xi;
}

Real Fluid3D::CipCsl3Z(Int32 i, Int32 j, Int32 k, Real v,
                       const Vector3DR &sia, const Vector3DR &via) noexcept
{
  auto xi = -v * dt_;

  Real c[4];
  Int32 kp, kv;
  Real dz;
  if(v < 0.0){
    kp = k + 1;
    kv = k;
    dz = z_[kp] - z_[k];
  }else{
    kp = k - 1;
    kv = k - 1;
    dz = z_[kp] - z_[k];
  }
  assert(std::fabs(xi) <= std::fabs(dz));
  auto d = GradientZ(z_, sia, via, i, j, kv);
  c[0] = sia(i, j, k);
  c[1] = 2.0 / dz * (3.0 * via(i, j, kv) - 3.0 * sia(i, j, k) - d * dz);
  c[2] = 3.0 / (dz * dz)
       * (-2.0 * via(i, j, kv) + 3.0 * sia(i, j, k) - sia(i, j, kp) + 2.0 * d * dz);
  c[3] = 4.0 / (dz * dz * dz) * (sia(i, j, kp) - sia(i, j, k) - d * dz);

  return c[0] + c[1] * xi + c[2] * xi * xi + c[3] * xi * xi * xi;
}
