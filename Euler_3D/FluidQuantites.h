#ifndef FLUIDQUANTITES_H
#define FLUIDQUANTITES_H

#include <type_traits>

#include "global_definitions.h"
#include "Vector3D.h"
class FluidQuantity3D
{
public:
  FluidQuantity3D() = default;
  FluidQuantity3D(Real d, Real u, Real v, Real w, Real p) noexcept
    : d_(d), u_(u), v_(v), w_(w), p_(p) {}

  FluidQuantity3D(const FluidQuantity3D&) = default;
  FluidQuantity3D& operator=(const FluidQuantity3D&) = default;

  FluidQuantity3D(FluidQuantity3D&&) = default;
  FluidQuantity3D& operator=(FluidQuantity3D&&) = default;

  ~FluidQuantity3D() = default;

  const Real& d() const noexcept{ return d_; }
  const Real& u() const noexcept{ return u_; }
  const Real& v() const noexcept{ return v_; }
  const Real& w() const noexcept{ return w_; }
  const Real& p() const noexcept{ return p_; }
  Real& d() noexcept{ return d_; }
  Real& u() noexcept{ return u_; }
  Real& v() noexcept{ return v_; }
  Real& w() noexcept{ return w_; }
  Real& p() noexcept{ return p_; }
  void d(Real val) noexcept{ d_ = val; }
  void u(Real val) noexcept{ u_ = val; }
  void v(Real val) noexcept{ v_ = val; }
  void w(Real val) noexcept{ w_ = val; }
  void p(Real val) noexcept{ p_ = val; }

  FluidQuantity3D& operator+=(const FluidQuantity3D& other)
  {
    d_ += other.d_;
    u_ += other.u_;
    v_ += other.v_;
    w_ += other.w_;
    p_ += other.p_;
    return *this;
  }

  FluidQuantity3D& operator-=(const FluidQuantity3D& other)
  {
    d_ -= other.d_;
    u_ -= other.u_;
    v_ -= other.v_;
    w_ -= other.w_;
    p_ -= other.p_;
    return *this;
  }

  FluidQuantity3D& operator*=(Real x)
  {
    d_ *= x;
    u_ *= x;
    v_ *= x;
    w_ *= x;
    p_ *= x;
    return *this;
  }
private:
  Real d_;
  Real u_;
  Real v_;
  Real w_;
  Real p_;
};

inline FluidQuantity3D operator+(const FluidQuantity3D &x, const FluidQuantity3D &y)
{
  FluidQuantity3D ret(x);
  ret += y;
  return ret;
}

inline FluidQuantity3D operator-(const FluidQuantity3D &x, const FluidQuantity3D &y)
{
  FluidQuantity3D ret(x);
  ret -= y;
  return ret;
}

inline FluidQuantity3D operator*(const FluidQuantity3D &x, Real y)
{
  FluidQuantity3D ret(x);
  return ret *= y;
}

inline FluidQuantity3D operator*(Real x, const FluidQuantity3D &y)
{
  return y * x;
}

class FluidQuantity3DArray
{
public:
  FluidQuantity3DArray() = default;
  FluidQuantity3DArray(Int32 nx, Int32 ny, Int32 nz)
    : data_(nx, ny, nz) {}

  FluidQuantity3DArray(const FluidQuantity3DArray&) = delete;
  FluidQuantity3DArray& operator=(const FluidQuantity3DArray&) = delete;

  FluidQuantity3DArray(FluidQuantity3DArray&&) = default;
  FluidQuantity3DArray& operator=(FluidQuantity3DArray&&) = default;

  ~FluidQuantity3DArray() = default;

  const FluidQuantity3D& operator()(Int32 i, Int32 j, Int32 k) const noexcept
  { return data_(i, j, k); }
  FluidQuantity3D& operator()(Int32 i, Int32 j, Int32 k) noexcept
  { return data_(i, j, k); }
  void operator ()(const FluidQuantity3D& val, Int32 i, Int32 j, Int32 k) noexcept
  { data_(i, j, k) = val; }

  const Real& d(Int32 i, Int32 j, Int32 k) const noexcept
  { return data_(i, j, k).d(); }
  const Real& u(Int32 i, Int32 j, Int32 k) const noexcept
  { return data_(i, j, k).u(); }
  const Real& v(Int32 i, Int32 j, Int32 k) const noexcept
  { return data_(i, j, k).v(); }
  const Real& w(Int32 i, Int32 j, Int32 k) const noexcept
  { return data_(i, j, k).w(); }
  const Real& p(Int32 i, Int32 j, Int32 k) const noexcept
  { return data_(i, j, k).p(); }

  Real& d(Int32 i, Int32 j, Int32 k) noexcept
  { return data_(i, j, k).d(); }
  Real& u(Int32 i, Int32 j, Int32 k) noexcept
  { return data_(i, j, k).u(); }
  Real& v(Int32 i, Int32 j, Int32 k) noexcept
  { return data_(i, j, k).v(); }
  Real& w(Int32 i, Int32 j, Int32 k) noexcept
  { return data_(i, j, k).w(); }
  Real& p(Int32 i, Int32 j, Int32 k) noexcept
  { return data_(i, j, k).p(); }

  void d(Real val, Int32 i, Int32 j, Int32 k) noexcept
  { data_(i, j, k).d(val); }
  void u(Real val, Int32 i, Int32 j, Int32 k) noexcept
  { data_(i, j, k).u(val); }
  void v(Real val, Int32 i, Int32 j, Int32 k) noexcept
  { data_(i, j, k).v(val); }
  void w(Real val, Int32 i, Int32 j, Int32 k) noexcept
  { data_(i, j, k).w(val); }
  void p(Real val, Int32 i, Int32 j, Int32 k) noexcept
  { data_(i, j, k).p(val); }

  Int32 size_x() const noexcept
  { return data_.size_x(); }
  Int32 size_y() const noexcept
  { return data_.size_y(); }
  Int32 size_z() const noexcept
  { return data_.size_z(); }
  Int32 size() const noexcept
  { return data_.size(); }

  void swap(FluidQuantity3DArray& other) noexcept
  { data_.swap(other.data_); }

private:
  Vector3D<FluidQuantity3D> data_;
};

inline void swap(FluidQuantity3DArray &a, FluidQuantity3DArray &b) noexcept
{
  a.swap(b);
}

struct FQAccessTagBase {};

struct TagD : public FQAccessTagBase {};
struct TagU : public FQAccessTagBase {};
struct TagV : public FQAccessTagBase {};
struct TagW : public FQAccessTagBase {};
struct TagP : public FQAccessTagBase {};

template <typename Tag>
class FQAccessor
{
public:
  explicit FQAccessor(const FluidQuantity3DArray &v)
    : value_(v)
  {
    static_assert(std::is_base_of<FQAccessTagBase, Tag>::value,
                  "Tag is derived from FQAccessTagBase!!");
  }
  FQAccessor() = delete;
  FQAccessor(const FQAccessor&) = default;
  FQAccessor& operator=(const FQAccessor&) = default;
  FQAccessor(FQAccessor&&) = default;
  FQAccessor& operator =(FQAccessor&&) = default;
  ~FQAccessor() = default;

  const Real & operator()(Int32 i, Int32 j, Int32 k) const noexcept
  {
    return this->operator ()(i, j, k, tag_);
  }
private:
  Tag tag_;
  const FluidQuantity3DArray &value_;

  const Real& operator()(Int32 i, Int32 j, Int32 k, TagD) const noexcept
  {
    return value_.d(i, j, k);
  }
  const Real& operator()(Int32 i, Int32 j, Int32 k, TagU) const noexcept
  {
    return value_.u(i, j, k);
  }
  const Real& operator()(Int32 i, Int32 j, Int32 k, TagV) const noexcept
  {
    return value_.v(i, j, k);
  }
  const Real& operator()(Int32 i, Int32 j, Int32 k, TagW) const noexcept
  {
    return value_.w(i, j, k);
  }
  const Real& operator()(Int32 i, Int32 j, Int32 k, TagP) const noexcept
  {
    return value_.p(i, j, k);
  }
};

#endif // FLUIDQUANTITES_H
