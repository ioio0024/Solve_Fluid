#include <cassert>
#include <cmath>

#include "global_definitions.h"
#include "FluidQuantites.h"

inline Real Sign(Real x) {
  if(x < 0.0){
    return -1.0;
  } else {
    return 1.0;
  }
}

inline Real Min(Real a, Real b, Real c){
  return std::fmin(std::fmin(a, b), c);
}

class CipCsl3X
{
public:
  template <typename Tag>
  Real operator()(Int32 i, Int32 j, Int32 k, Real xi, const Vector3DR &x,
                  FQAccessor<Tag> sia, FQAccessor<Tag> via) noexcept;
private:
  template <typename Tag>
  Real Gradient(const Vector3DR &x, FQAccessor<Tag> sia, FQAccessor<Tag> via,
                Int32 i, Int32 j, Int32 k) noexcept;
};

template <typename Tag>
Real CipCsl3X::operator()(Int32 i, Int32 j, Int32 k, Real xi, const Vector3DR &x,
                          FQAccessor<Tag> sia, FQAccessor<Tag> via) noexcept
{
  Real c[4];
  Int32 ip, iv;
  if(xi >= 0.0){
    ip = i + 1;
    iv = i;
  }else{
    ip = i - 1;
    iv = i - 1;
  }
  Real dx = x(ip) - x(i);
  assert(std::fabs(xi) <= std::fabs(dx));
  auto d = Gradient(x, sia, via, iv, j, k);
  c[0] = sia(i, j, k);
  c[1] = 2.0 / dx * (3.0 * via(iv, j, k) - 3.0 * sia(i, j, k) - d * dx);
  c[2] = 3.0 / (dx * dx)
       * (-2.0 * via(iv, j, k) + 3.0 * sia(i, j, k) - sia(ip, j, k) + 2.0 * d * dx);
  c[3] = 4.0 / (dx * dx * dx) * (sia(ip, j, k) - sia(i, j, k) - d * dx);

  return c[0] + c[1] * xi + c[2] * xi * xi + c[3] * xi * xi * xi;
}

template <typename Tag>
Real CipCsl3X::Gradient(const Vector3DR &x, FQAccessor<Tag> sia, FQAccessor<Tag> via,
                        Int32 i, Int32 j, Int32 k) noexcept
{
  Real dx[3];
  Real f[3];

  for(auto n = 0; n < 3; ++n){
    auto iv = i - 1+ n;
    dx[n] = x(iv + 1) - x(iv);
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


class CipCsl3Y
{
public:
  template <typename Tag>
  Real operator()(Int32 i, Int32 j, Int32 k, Real xi, const Vector3DR &x,
                  FQAccessor<Tag> sia, FQAccessor<Tag> via) noexcept;
private:
  template <typename Tag>
  Real Gradient(const Vector3DR &x, FQAccessor<Tag> sia, FQAccessor<Tag> via,
                Int32 i, Int32 j, Int32 k) noexcept;
};


template <typename Tag>
Real CipCsl3Y::operator()(Int32 i, Int32 j, Int32 k, Real xi, const Vector3DR &y,
                          FQAccessor<Tag> sia, FQAccessor<Tag> via) noexcept
{
  Real c[4];
  Int32 jp, jv;
  if(xi >= 0.0){
    jp = j + 1;
    jv = j;
  }else{
    jp = j - 1;
    jv = j - 1;
  }
  Real dy = y(jp) - y(j);
  assert(std::fabs(xi) <= std::fabs(dy));
  auto d = Gradient(y, sia, via, i, jv, k);
  c[0] = sia(i, j, k);
  c[1] = 2.0 / dy * (3.0 * via(i, jv, k) - 3.0 * sia(i, j, k) - d * dy);
  c[2] = 3.0 / (dy * dy)
       * (-2.0 * via(i, jv, k) + 3.0 * sia(i, j, k) - sia(i, jp, k) + 2.0 * d * dy);
  c[3] = 4.0 / (dy * dy * dy) * (sia(i, jp, k) - sia(i, j, k) - d * dy);

  return c[0] + c[1] * xi + c[2] * xi * xi + c[3] * xi * xi * xi;
}

template <typename Tag>
Real CipCsl3Y::Gradient(const Vector3DR &x, FQAccessor<Tag> sia, FQAccessor<Tag> via,
                        Int32 i, Int32 j, Int32 k) noexcept
{
  Real dx[3];
  Real f[3];

  for(auto n = 0; n < 3; ++n){
    auto jv = j - 1+ n;
    dx[n] = x(jv + 1) - x(jv);
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


class CipCsl3Z
{
public:
  template <typename Tag>
  Real operator()(Int32 i, Int32 j, Int32 k, Real xi, const Vector3DR &x,
                  FQAccessor<Tag> sia, FQAccessor<Tag> via) noexcept;
private:
  template <typename Tag>
  Real Gradient(const Vector3DR &x, FQAccessor<Tag> sia, FQAccessor<Tag> via,
                Int32 i, Int32 j, Int32 k) noexcept;
};


template <typename Tag>
Real CipCsl3Z::operator()(Int32 i, Int32 j, Int32 k, Real xi, const Vector3DR &z,
                          FQAccessor<Tag> sia, FQAccessor<Tag> via) noexcept
{
  Real c[4];
  Int32 kp, kv;
  if(xi >= 0.0){
    kp = k + 1;
    kv = k;
  }else{
    kp = k - 1;
    kv = k - 1;
  }
  Real dz = z(kp) - z(k);
  assert(std::fabs(xi) <= std::fabs(dz));
  auto d = Gradient(z, sia, via, i, j, kv);
  c[0] = sia(i, j, k);
  c[1] = 2.0 / dz * (3.0 * via(i, j, kv) - 3.0 * sia(i, j, k) - d * dz);
  c[2] = 3.0 / (dz * dz)
       * (-2.0 * via(i, j, kv) + 3.0 * sia(i, j, k) - sia(i, j, kp) + 2.0 * d * dz);
  c[3] = 4.0 / (dz * dz * dz) * (sia(i, j, kp) - sia(i, j, k) - d * dz);

  return c[0] + c[1] * xi + c[2] * xi * xi + c[3] * xi * xi * xi;
}

template <typename Tag>
Real CipCsl3Z::Gradient(const Vector3DR &x, FQAccessor<Tag> sia, FQAccessor<Tag> via,
                        Int32 i, Int32 j, Int32 k) noexcept
{
  Real dx[3];
  Real f[3];

  for(auto n = 0; n < 3; ++n){
    auto kv = k - 1+ n;
    dx[n] = x(kv + 1) - x(kv);
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
