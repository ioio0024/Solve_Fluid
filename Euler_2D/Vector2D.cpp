#include <utility>

#include "Vector2D.h"

Vector2D::Vector2D(Int32 size_x, Int32 size_y)
  : size_x_{size_x}, size_y_{size_y},
    size_{size_x_ * size_y_}, data_(size_) {}

Vector2D::Vector2D(Int32 size_x, Int32 size_y, Real val)
  : size_x_{size_x}, size_y_{size_y},
    size_{size_x_ * size_y_}, data_(size_, val) {}

void Vector2D::resize(Int32 size_x, Int32 size_y)
{
  size_x_ = size_x;
  size_y_ = size_y;
  size_ = size_x_ * size_y_;

  data_.resize(size_);
}

void Vector2D::swap(Vector2D &x) noexcept
{
  std::swap(size_x_, x.size_x_);
  std::swap(size_y_, x.size_y_);
  std::swap(size_, x.size_);
  std::swap(data_, x.data_);
}

void swap(Vector2D &x, Vector2D &y) noexcept
{
  x.swap(y);
}
