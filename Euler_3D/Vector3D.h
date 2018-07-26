#ifndef VECTOR3D_H
#define VECTOR3D_H

#include "global_definitions.h"

#include <cassert>
#include <utility>

template <typename T>
class Vector3D
{
public:
  Vector3D();
  Vector3D(Int32 nx, Int32 ny = 1, Int32 nz = 1); // 0で初期化
  Vector3D(T val, Int32 nx, Int32 ny = 1, Int32 nz = 1); // valで初期化

  ~Vector3D(){
    delete[] data_;
  }

  Vector3D(const Vector3D&) = delete;
  Vector3D& operator= (const Vector3D&) = delete;

  Vector3D(Vector3D&& other);
  Vector3D& operator= (Vector3D&& other);

  Int32 size_x() const noexcept{
    return nx_;
  }
  Int32 size_y() const noexcept{
    return ny_;
  }
  Int32 size_z() const noexcept{
    return nz_;
  }
  Int32 size() const noexcept{
    return nx_ * ny_ * nz_;
  }

  const T& operator()(Int32 i, Int32 j = 0, Int32 k = 0) const noexcept
  {
    assert(i >= 0 && i < nx_);
    assert(j >= 0 && j < ny_);
    assert(k >= 0 && k < nz_);
    return data_[nx_ * (ny_ * k + j) + i];
  }

  T& operator()(Int32 i, Int32 j = 0, Int32 k = 0) noexcept
  {
    assert(i >= 0 && i < nx_);
    assert(j >= 0 && j < ny_);
    assert(k >= 0 && k < nz_);
    return data_[nx_ * (ny_ * k + j) + i];
  }


  void swap(Vector3D& other);

private:
  Int32 nx_;
  Int32 ny_;
  Int32 nz_;
  T *data_;
};

template <typename T>
Vector3D<T>::Vector3D()
  : nx_(1), ny_(1), nz_(1), data_(nullptr)
{}

template <typename T>
Vector3D<T>::Vector3D(Int32 nx, Int32 ny, Int32 nz)
 : nx_(nx), ny_(ny), nz_(nz), data_(new T[nx * ny * nz]())
{}

template <typename T>
Vector3D<T>::Vector3D(T val, Int32 nx, Int32 ny, Int32 nz)
 : nx_(nx), ny_(ny), nz_(nz), data_(new T[nx * ny * nz])
{
  for(auto n = 0; n < nx * ny * nz; ++n){
    data_[n] = val;
  }
}

template <typename T>
Vector3D<T>::Vector3D(Vector3D &&other)
 : nx_(other.nx_), ny_(other.ny_), nz_(other.nz_), data_(nullptr)
{
  std::swap(data_, other.data_);
}

template <typename T>
Vector3D<T>& Vector3D<T>::operator =(Vector3D &&other)
{
  nx_ = other.nx_;
  ny_ = other.ny_;
  nz_ = other.nz_;
  std::swap(data_, other.data_);
  return *this;
}

template <typename T>
void Vector3D<T>::swap(Vector3D &other)
{
  std::swap(nx_, other.nx_);
  std::swap(ny_, other.ny_);
  std::swap(nz_, other.nz_);
  std::swap(data_, other.data_);
}

template <typename T>
void swap(Vector3D<T> &l, Vector3D<T> &r)
{
  l.swap(r);
}

//extern template class Vector3D<Real>;
//extern template class Vector3D<int>;

using Vector3DR = Vector3D<Real>;
using Vector3DI = Vector3D<int>;

#endif // VECTOR3D_H
