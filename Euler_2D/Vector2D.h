#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <cassert>
#include <vector>

#include "global_definitions.h"

class Vector2D
{
public:
  Vector2D() = default;
  Vector2D(Int32 size_x, Int32 size_y);
  Vector2D(Int32 size_x, Int32 size_y, Real val);

  Vector2D(const Vector2D &) = default;
  Vector2D& operator =(const Vector2D &) = default;

  Vector2D(Vector2D&&) = default;
  Vector2D& operator =(Vector2D&&) = default;

  ~Vector2D() = default;

  Int32 size_x() const noexcept { return size_x_; }
  Int32 size_y() const noexcept { return size_y_; }
  Int32 size() const noexcept { return size_; }

  const Real& operator ()(Int32 i, Int32 j) const noexcept
  {
    assert(i >= 0 && i < size_x_);
    assert(j >= 0 && j < size_y_);
    return data_[j * size_x_ + i];
  }
  Real& operator ()(Int32 i, Int32 j) noexcept
  {
    assert(i >= 0 && i < size_x_);
    assert(j >= 0 && j < size_y_);
    return data_[j * size_x_ + i];
  }

  void resize(Int32 size_x, Int32 size_y);

  void swap(Vector2D &x) noexcept;

private:
  Int32 size_x_;
  Int32 size_y_;
  Int32 size_;
  Rvec data_;
};

void swap(Vector2D &x, Vector2D &y) noexcept;

#endif // VECTOR2D_H
