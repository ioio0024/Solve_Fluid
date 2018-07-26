template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
                ::Initialize(Int32 size_x, Int32 size_y, Int32 size_z,
                             Real xmin, Real xmax,
                             Real ymin, Real ymax,
                             Real zmin, Real zmax,
                             Real gamma,
                             std::function<void (const Vector3DR &,
                                                 const Vector3DR &,
                                                 const Vector3DR &,
                                                 FluidQuantity3DArray *
                                                 )> init_func)
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

  x_ = Vector3DR(nx_ + 1);
  y_ = Vector3DR(ny_ + 1);
  z_ = Vector3DR(nz_ + 1);

  dt_ = 1.0e3;
  time_ = 0.0;
  gamma_ = gamma;

  pv_ = FluidQuantity3DArray(nx_ + 1, ny_ + 1, nz_ + 1);

  lia_x_ = FluidQuantity3DArray(nx_    , ny_ + 1, nz_ + 1);
  lia_y_ = FluidQuantity3DArray(nx_ + 1, ny_    , nz_ + 1);
  lia_z_ = FluidQuantity3DArray(nx_ + 1, ny_ + 1, nz_    );

  sia_x_ = FluidQuantity3DArray(nx_ + 1, ny_    , nz_    );
  sia_y_ = FluidQuantity3DArray(nx_    , ny_ + 1, nz_    );
  sia_z_ = FluidQuantity3DArray(nx_    , ny_    , nz_ + 1);

  via_ = FluidQuantity3DArray(nx_, ny_, nz_);

  pv_new_ = FluidQuantity3DArray(nx_ + 1, ny_ + 1, nz_ + 1);

  lia_x_new_ = FluidQuantity3DArray(nx_    , ny_ + 1, nz_ + 1);
  lia_y_new_ = FluidQuantity3DArray(nx_ + 1, ny_    , nz_ + 1);
  lia_z_new_ = FluidQuantity3DArray(nx_ + 1, ny_ + 1, nz_    );

  sia_x_new_ = FluidQuantity3DArray(nx_ + 1, ny_    , nz_    );
  sia_y_new_ = FluidQuantity3DArray(nx_    , ny_ + 1, nz_    );
  sia_z_new_ = FluidQuantity3DArray(nx_    , ny_    , nz_ + 1);

  via_new_ = FluidQuantity3DArray(nx_, ny_, nz_);

  flux_ = FluidQuantity3DArray(nx_ + 1, ny_ + 1, nz_ + 1);

  auto dx = (xmax_ - xmin_) / size_x;
  auto dy = (ymax_ - ymin_) / size_y;
  auto dz = (zmax_ - zmin_) / size_z;

  for(auto i = 0; i < nx_ + 1; ++i){
    x_(i) = xmin_ + dx * (i - GHOST);
  }
  for(auto j = 0; j < ny_ + 1; ++j){
    y_(j) = ymin_ + dy * (j - GHOST);
  }
  for(auto k = 0; k < nz_ + 1; ++k){
    z_(k) = zmin_ + dz * (k - GHOST);
  }

  init_func(x_, y_, z_, &pv_);

  // LIA x
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        lia_x_(i, j, k) = 0.5 * (pv_(i + 1, j, k) + pv_(i, j, k));
      }
    }
  }

  // LIA y
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        lia_y_(i, j, k) = 0.5 * (pv_(i, j + 1, k) + pv_(i, j, k));
      }
    }
  }

  // LIA z
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        lia_z_(i, j, k) = 0.5 * (pv_(i, j, k + 1) + pv_(i, j, k));
      }
    }
  }

  // SIA x
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        sia_x_(i, j, k) = 0.25 * (pv_(i, j + 1, k + 1) + pv_(i, j + 1, k)
                                + pv_(i, j, k + 1) + pv_(i, j, k));
      }
    }
  }

  // SIA y
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        sia_y_(i, j, k) = 0.25 * (pv_(i + 1, j, k + 1) + pv_(i + 1, j, k)
                                 + pv_(i, j, k + 1) + pv_(i, j, k));
      }
    }
  }

  // SIA z
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_; ++i){
        sia_z_(i, j, k) = 0.25 * (pv_(i + 1, j + 1, k) + pv_(i + 1, j, k)
                                  + pv_(i, j + 1, k) + pv_(i, j, k));
      }
    }
  }

  // VIA
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_; ++i){
        via_(i, j, k) = 0.125 * (pv_(i + 1, j + 1, k + 1) + pv_(i + 1, j + 1, k)
                                 + pv_(i, j + 1, k + 1) + pv_(i + 1, j, k + 1)
                                 + pv_(i + 1, j, k) + pv_(i, j + 1, k)
                                 + pv_(i, j, k + 1) + pv_(i, j, k));
      }
    }
  }

  NewDt();
  dt_ *= 0.01;
}
