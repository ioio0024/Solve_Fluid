#include "global_definitions.h"
#include "Fluid3D.h"

Fluid3D& Fluid3D::Initialize(Int32 size_x, Int32 size_y, Int32 size_z,
                             Real xmin, Real xmax,
                             Real ymin, Real ymax,
                             Real zmin, Real zmax, Real gamma,
                             std::function<void (const Rvec &,
                                                 const Rvec &,
                                                 const Rvec &,
                                                 Vector3DR *,
                                                 Vector3DR *,
                                                 Vector3DR *,
                                                 Vector3DR *,
                                                 Vector3DR *)>
                             init_func)
{
  nx_ = size_x + 2 * GHOST;
  ny_ = size_y + 2 * GHOST;
  nz_ = size_z + 2 * GHOST;
  xmin_ = xmin;
  xmax_ = xmax;
  ymin_ = ymin;
  ymax_ = ymax;
  zmin_ = zmin;
  zmax_ = zmax;

  x_.resize(nx_ + 1);
  y_.resize(ny_ + 1);
  z_.resize(nz_ + 1);

  dt_ = 1.0e3;
  time_ = 0.0;
  gamma_ = gamma;

  density_pv_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  vx_pv_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  vy_pv_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  vz_pv_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  pressure_pv_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);

  density_lia_x_ = Vector3DR(nx_, ny_ + 1, nz_ + 1);
  vx_lia_x_ = Vector3DR(nx_, ny_ + 1, nz_ + 1);
  vy_lia_x_ = Vector3DR(nx_, ny_ + 1, nz_ + 1);
  vz_lia_x_ = Vector3DR(nx_, ny_ + 1, nz_ + 1);
  pressure_lia_x_ = Vector3DR(nx_, ny_ + 1, nz_ + 1);

  density_lia_y_ = Vector3DR(nx_ + 1, ny_, nz_ + 1);
  vx_lia_y_ = Vector3DR(nx_ + 1, ny_, nz_ + 1);
  vy_lia_y_ = Vector3DR(nx_ + 1, ny_, nz_ + 1);
  vz_lia_y_ = Vector3DR(nx_ + 1, ny_, nz_ + 1);
  pressure_lia_y_ = Vector3DR(nx_ + 1, ny_, nz_ + 1);

  density_lia_z_ = Vector3DR(nx_ + 1, ny_ + 1, nz_);
  vx_lia_z_ = Vector3DR(nx_ + 1, ny_ + 1, nz_);
  vy_lia_z_ = Vector3DR(nx_ + 1, ny_ + 1, nz_);
  vz_lia_z_ = Vector3DR(nx_ + 1, ny_ + 1, nz_);
  pressure_lia_z_ = Vector3DR(nx_ + 1, ny_ + 1, nz_);

  density_sia_x_ = Vector3DR(nx_ + 1, ny_, nz_);
  vx_sia_x_ = Vector3DR(nx_ + 1, ny_, nz_);
  vy_sia_x_ = Vector3DR(nx_ + 1, ny_, nz_);
  vz_sia_x_ = Vector3DR(nx_ + 1, ny_, nz_);
  pressure_sia_x_ = Vector3DR(nx_ + 1, ny_, nz_);

  density_sia_y_ = Vector3DR(nx_, ny_ + 1, nz_);
  vx_sia_y_ = Vector3DR(nx_, ny_ + 1, nz_);
  vy_sia_y_ = Vector3DR(nx_, ny_ + 1, nz_);
  vz_sia_y_ = Vector3DR(nx_, ny_ + 1, nz_);
  pressure_sia_y_ = Vector3DR(nx_, ny_ + 1, nz_);

  density_sia_z_ = Vector3DR(nx_, ny_, nz_ + 1);
  vx_sia_z_ = Vector3DR(nx_, ny_, nz_ + 1);
  vy_sia_z_ = Vector3DR(nx_, ny_, nz_ + 1);
  vz_sia_z_ = Vector3DR(nx_, ny_, nz_ + 1);
  pressure_sia_z_ = Vector3DR(nx_, ny_, nz_ + 1);

  density_via_ = Vector3DR(nx_, ny_, nz_);
  vx_via_ = Vector3DR(nx_, ny_, nz_);
  vy_via_ = Vector3DR(nx_, ny_, nz_);
  vz_via_ = Vector3DR(nx_, ny_, nz_);
  pressure_via_ = Vector3DR(nx_, ny_, nz_);

  density_pv_new_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  vx_pv_new_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  vy_pv_new_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  vz_pv_new_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  pressure_pv_new_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);

  density_lia_x_new_ = Vector3DR(nx_, ny_ + 1, nz_ + 1);
  vx_lia_x_new_ = Vector3DR(nx_, ny_ + 1, nz_ + 1);
  vy_lia_x_new_ = Vector3DR(nx_, ny_ + 1, nz_ + 1);
  vz_lia_x_new_ = Vector3DR(nx_, ny_ + 1, nz_ + 1);
  pressure_lia_x_new_ = Vector3DR(nx_, ny_ + 1, nz_ + 1);

  density_lia_y_new_ = Vector3DR(nx_ + 1, ny_, nz_ + 1);
  vx_lia_y_new_ = Vector3DR(nx_ + 1, ny_, nz_ + 1);
  vy_lia_y_new_ = Vector3DR(nx_ + 1, ny_, nz_ + 1);
  vz_lia_y_new_ = Vector3DR(nx_ + 1, ny_, nz_ + 1);
  pressure_lia_y_new_ = Vector3DR(nx_ + 1, ny_, nz_ + 1);

  density_lia_z_new_ = Vector3DR(nx_ + 1, ny_ + 1, nz_);
  vx_lia_z_new_ = Vector3DR(nx_ + 1, ny_ + 1, nz_);
  vy_lia_z_new_ = Vector3DR(nx_ + 1, ny_ + 1, nz_);
  vz_lia_z_new_ = Vector3DR(nx_ + 1, ny_ + 1, nz_);
  pressure_lia_z_new_ = Vector3DR(nx_ + 1, ny_ + 1, nz_);

  density_sia_x_new_ = Vector3DR(nx_ + 1, ny_, nz_);
  vx_sia_x_new_ = Vector3DR(nx_ + 1, ny_, nz_);
  vy_sia_x_new_ = Vector3DR(nx_ + 1, ny_, nz_);
  vz_sia_x_new_ = Vector3DR(nx_ + 1, ny_, nz_);
  pressure_sia_x_new_ = Vector3DR(nx_ + 1, ny_, nz_);

  density_sia_y_new_ = Vector3DR(nx_, ny_ + 1, nz_);
  vx_sia_y_new_ = Vector3DR(nx_, ny_ + 1, nz_);
  vy_sia_y_new_ = Vector3DR(nx_, ny_ + 1, nz_);
  vz_sia_y_new_ = Vector3DR(nx_, ny_ + 1, nz_);
  pressure_sia_y_new_ = Vector3DR(nx_, ny_ + 1, nz_);

  density_sia_z_new_ = Vector3DR(nx_, ny_, nz_ + 1);
  vx_sia_z_new_ = Vector3DR(nx_, ny_, nz_ + 1);
  vy_sia_z_new_ = Vector3DR(nx_, ny_, nz_ + 1);
  vz_sia_z_new_ = Vector3DR(nx_, ny_, nz_ + 1);
  pressure_sia_z_new_ = Vector3DR(nx_, ny_, nz_ + 1);

  density_via_new_ = Vector3DR(nx_, ny_, nz_);
  vx_via_new_ = Vector3DR(nx_, ny_, nz_);
  vy_via_new_ = Vector3DR(nx_, ny_, nz_);
  vz_via_new_ = Vector3DR(nx_, ny_, nz_);
  pressure_via_new_ = Vector3DR(nx_, ny_, nz_);

  dflux_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  lflux_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  mflux_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  nflux_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);
  eflux_ = Vector3DR(nx_ + 1, ny_ + 1, nz_ + 1);

  auto dx = (xmax_ - xmin_) / size_x;
  auto dy = (ymax_ - ymin_) / size_y;
  auto dz = (zmax_ - zmin_) / size_z;

  for(auto i = 0; i < nx_ + 1; ++i){
    x_[i] = xmin_ + dx * (i - GHOST);
  }
  for(auto j = 0; j < ny_ + 1; ++j){
    y_[j] = ymin_ + dy * (j - GHOST);
  }
  for(auto k = 0; k < nz_ + 1; ++k){
    z_[k] = zmin_ + dz * (k - GHOST);
  }

  init_func(x_, y_, z_,
            &density_pv_,
            &vx_pv_,
            &vy_pv_,
            &vz_pv_,
            &pressure_pv_);

  // LIA x
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_lia_x_(i, j, k) = 0.5 * (density_pv_(i + 1,j, k)
                                         + density_pv_(i, j, k));
        vx_lia_x_(i, j, k) = 0.5 * (vx_pv_(i + 1, j, k)
                                    + vx_pv_(i, j, k));

        vy_lia_x_(i, j, k) = 0.5 * (vy_pv_(i + 1, j, k)
                                    + vy_pv_(i, j, k));
        vz_lia_x_(i, j, k) = 0.5 * (vz_pv_(i + 1, j, k)
                                    + vz_pv_(i, j, k));
        pressure_lia_x_(i, j, k) = 0.5 * (pressure_pv_(i + 1, j, k)
                                          + pressure_pv_(i, j, k));
      }
    }
  }

  // LIA y
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        density_lia_y_(i, j, k) = 0.5 * (density_pv_(i, j + 1, k)
                                         + density_pv_(i, j, k));
        vx_lia_y_(i, j, k) = 0.5 * (vx_pv_(i, j + 1, k)
                                    + vx_pv_(i, j, k));

        vy_lia_y_(i, j, k) = 0.5 * (vy_pv_(i, j + 1, k)
                                    + vy_pv_(i, j, k));
        vz_lia_y_(i, j, k) = 0.5 * (vz_pv_(i, j + 1, k)
                                    + vz_pv_(i, j, k));
        pressure_lia_y_(i, j, k) = 0.5 * (pressure_pv_(i, j + 1, k)
                                          + pressure_pv_(i, j, k));
      }
    }
  }

  // LIA z
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        density_lia_z_(i, j, k) = 0.5 * (density_pv_(i, j, k + 1)
                                         + density_pv_(i, j, k));
        vx_lia_z_(i, j, k) = 0.5 * (vx_pv_(i, j, k + 1)
                                    + vx_pv_(i, j, k));

        vy_lia_z_(i, j, k) = 0.5 * (vy_pv_(i, j, k + 1)
                                    + vy_pv_(i, j, k));
        vz_lia_z_(i, j, k) = 0.5 * (vz_pv_(i, j, k + 1)
                                    + vz_pv_(i, j, k));
        pressure_lia_z_(i, j, k) = 0.5 * (pressure_pv_(i, j, k + 1)
                                          + pressure_pv_(i, j, k));
      }
    }
  }

  // SIA x
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        density_sia_x_(i, j, k) = 0.25 * (density_pv_(i,j + 1, k + 1)
                                         + density_pv_(i, j + 1, k)
                                         + density_pv_(i, j, k + 1)
                                         + density_pv_(i, j, k));
        vx_sia_x_(i, j, k) = 0.25 * (vx_pv_(i, j + 1, k + 1)
                                    + vx_pv_(i, j + 1, k)
                                    + vx_pv_(i, j, k + 1)
                                    + vx_pv_(i, j, k));
        vy_sia_x_(i, j, k) = 0.25 * (vy_pv_(i, j + 1, k + 1)
                                    + vy_pv_(i, j + 1, k)
                                    + vy_pv_(i, j, k + 1)
                                    + vy_pv_(i, j, k));
        vz_sia_x_(i, j, k) = 0.25 * (vz_pv_(i, j + 1, k + 1)
                                    + vz_pv_(i, j + 1, k)
                                    + vz_pv_(i, j, k + 1)
                                    + vz_pv_(i, j, k));
        pressure_sia_x_(i, j, k) = 0.25 * (pressure_pv_(i, j + 1, k + 1)
                                          + pressure_pv_(i, j + 1, k)
                                          + pressure_pv_(i, j, k + 1)
                                          + pressure_pv_(i, j, k));
      }
    }
  }

  // SIA y
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_sia_y_(i, j, k) = 0.25 * (density_pv_(i + 1,j, k + 1)
                                         + density_pv_(i + 1, j, k)
                                         + density_pv_(i, j, k + 1)
                                         + density_pv_(i, j, k));
        vx_sia_y_(i, j, k) = 0.25 * (vx_pv_(i + 1, j, k + 1)
                                    + vx_pv_(i + 1, j, k)
                                    + vx_pv_(i, j, k + 1)
                                    + vx_pv_(i, j, k));
        vy_sia_y_(i, j, k) = 0.25 * (vy_pv_(i + 1, j, k + 1)
                                    + vy_pv_(i + 1, j, k)
                                    + vy_pv_(i, j, k + 1)
                                    + vy_pv_(i, j, k));
        vz_sia_y_(i, j, k) = 0.25 * (vz_pv_(i + 1, j, k + 1)
                                    + vz_pv_(i + 1, j, k)
                                    + vz_pv_(i, j, k + 1)
                                    + vz_pv_(i, j, k));
        pressure_sia_y_(i, j, k) = 0.25 * (pressure_pv_(i + 1, j, k + 1)
                                          + pressure_pv_(i + 1, j, k)
                                          + pressure_pv_(i, j, k + 1)
                                          + pressure_pv_(i, j, k));
      }
    }
  }

  // SIA z
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_sia_z_(i, j, k) = 0.25 * (density_pv_(i + 1,j + 1, k)
                                         + density_pv_(i + 1, j, k)
                                         + density_pv_(i, j + 1, k)
                                         + density_pv_(i, j, k));
        vx_sia_z_(i, j, k) = 0.25 * (vx_pv_(i + 1, j + 1, k)
                                    + vx_pv_(i + 1, j, k)
                                    + vx_pv_(i, j + 1, k)
                                    + vx_pv_(i, j, k));
        vy_sia_z_(i, j, k) = 0.25 * (vy_pv_(i + 1, j + 1, k)
                                    + vy_pv_(i + 1, j, k)
                                    + vy_pv_(i, j + 1, k)
                                    + vy_pv_(i, j, k));
        vz_sia_z_(i, j, k) = 0.25 * (vz_pv_(i + 1, j + 1, k)
                                    + vz_pv_(i + 1, j, k)
                                    + vz_pv_(i, j + 1, k)
                                    + vz_pv_(i, j, k));
        pressure_sia_z_(i, j, k) = 0.25 * (pressure_pv_(i + 1, j + 1, k)
                                          + pressure_pv_(i + 1, j, k)
                                          + pressure_pv_(i, j + 1, k)
                                          + pressure_pv_(i, j, k));
      }
    }
  }

  // VIA
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_via_(i, j, k) = 0.125 * (density_pv_(i + 1,j + 1, k + 1)
                                         + density_pv_(i + 1, j + 1, k)
                                         + density_pv_(i, j + 1, k + 1)
                                         + density_pv_(i + 1, j, k + 1)
                                         + density_pv_(i + 1, j, k)
                                         + density_pv_(i, j + 1, k)
                                         + density_pv_(i, j, k + 1)
                                         + density_pv_(i, j, k));
        vx_via_(i, j, k) = 0.25 * (vx_pv_(i + 1, j + 1, k + 1)
                                   + vx_pv_(i + 1, j + 1, k)
                                   + vx_pv_(i, j + 1, k + 1)
                                   + vx_pv_(i + 1, j, k + 1)
                                   + vx_pv_(i + 1, j, k)
                                   + vx_pv_(i, j + 1, k)
                                   + vx_pv_(i, j, k + 1)
                                   + vx_pv_(i, j, k));
        vy_via_(i, j, k) = 0.25 * (vy_pv_(i + 1, j + 1, k + 1)
                                   + vy_pv_(i + 1, j + 1, k)
                                   + vy_pv_(i, j + 1, k + 1)
                                   + vy_pv_(i + 1, j, k + 1)
                                   + vy_pv_(i + 1, j, k)
                                   + vy_pv_(i, j + 1, k)
                                   + vy_pv_(i, j, k + 1)
                                   + vy_pv_(i, j, k));
        vz_via_(i, j, k) = 0.25 * (vz_pv_(i + 1, j + 1, k + 1)
                                   + vz_pv_(i + 1, j + 1, k)
                                   + vz_pv_(i, j + 1, k + 1)
                                   + vz_pv_(i + 1, j, k + 1)
                                   + vz_pv_(i + 1, j, k)
                                   + vz_pv_(i, j + 1, k)
                                   + vz_pv_(i, j, k + 1)
                                   + vz_pv_(i, j, k));
        pressure_via_(i, j, k) = 0.25 * (pressure_pv_(i + 1, j + 1, k + 1)
                                         + pressure_pv_(i + 1, j + 1, k)
                                         + pressure_pv_(i, j + 1, k + 1)
                                         + pressure_pv_(i + 1, j, k + 1)
                                         + pressure_pv_(i + 1, j, k)
                                         + pressure_pv_(i, j + 1, k)
                                         + pressure_pv_(i, j, k + 1)
                                         + pressure_pv_(i, j, k));
      }
    }
  }

  BoundaryConditions();

  NewDt();
  dt_ *= 0.01;

  return *this;
}
