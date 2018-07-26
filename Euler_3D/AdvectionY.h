template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
            ::AdvectionY(Int32 nz, Int32 nx,
                         const FluidQuantity3DArray &sia,
                         const FluidQuantity3DArray &via,
                         FluidQuantity3DArray *sia_new,
                         FluidQuantity3DArray *via_new
                         ) noexcept
{
  FluidQuantity3D rk[4], flux[3];
  Real engy[3];

  for(auto k = GHOST; k < nz - GHOST; ++k){
    for(auto j = GHOST; j < ny_ + 1 - GHOST; ++j){
      for(auto i = GHOST; i < nx - GHOST; ++i){
        rk[0].d() = sia.d(i, j, k);
        rk[0].u() = sia.u(i, j, k);
        rk[0].v() = sia.v(i, j, k);
        rk[0].w() = sia.w(i, j, k);
        rk[0].p() = sia.p(i, j, k);

        AdvCipRK3Y(0, i, j, k, rk, sia, via);
        AdvCipRK3Y(1, i, j, k, rk, sia, via);
        AdvCipRK3Y(2, i, j, k, rk, sia, via);

        for(auto l = 0; l < 3; ++l){
          engy[l] = rk[l].p() + 0.5 * (gamma_ - 1.0)
                  * rk[l].d() * (rk[l].u() * rk[l].u()
                               + rk[l].v() * rk[l].v()
                               + rk[l].w() * rk[l].w());
        }

        for(auto l = 0; l < 3; ++l){
          auto d = rk[l].d();
          auto u = rk[l].u();
          auto v = rk[l].v();
          auto w = rk[l].w();
          auto p = rk[l].p();

          flux[l].d() = d * v;
          flux[l].u() = d * u * v;
          flux[l].v() = d * v * v + p;
          flux[l].w() = d * v * w;
          flux[l].p() = v * (engy[l] + (gamma_ - 1.0) * p);
        }

        sia_new->d(i, j, k) = rk[3].d();
        sia_new->u(i, j, k) = rk[3].u();
        sia_new->v(i, j, k) = rk[3].v();
        sia_new->w(i, j, k) = rk[3].w();
        sia_new->p(i, j, k) = rk[3].p();

        flux_.d(i, j, k) = (flux[0].d() + flux[1].d() + 4.0 * flux[2].d()) / 6.0;
        flux_.u(i, j, k) = (flux[0].u() + flux[1].u() + 4.0 * flux[2].u()) / 6.0;
        flux_.v(i, j, k) = (flux[0].v() + flux[1].v() + 4.0 * flux[2].v()) / 6.0;
        flux_.w(i, j, k) = (flux[0].w() + flux[1].w() + 4.0 * flux[2].w()) / 6.0;
        flux_.p(i, j, k) = (flux[0].p() + flux[1].p() + 4.0 * flux[2].p()) / 6.0;
      }
    }
  }

  // VIA
  for(auto k = GHOST; k < nz - GHOST; ++k){
    for(auto j = GHOST; j < ny_ - GHOST; ++j){
      for(auto i = GHOST; i < nx - GHOST; ++i){
        auto dy = y_(j + 1) - y_(j);
        via_new->d(i, j, k) = via.d(i, j, k) - dt_ / dy
            * (flux_.d(i, j + 1, k) - flux_.d(i, j, k));
        via_new->u(i, j, k) = (via.d(i, j, k) * via.u(i, j, k) - dt_ / dy
                         * (flux_.u(i, j + 1, k) - flux_.u(i, j, k)))
                         / via_new->d(i, j, k);
        via_new->v(i, j, k) = (via.d(i, j, k) * via.v(i, j, k) - dt_ / dy
                         * (flux_.v(i, j + 1, k) - flux_.v(i, j, k)))
                         / via_new->d(i, j, k);
        via_new->w(i, j, k) = (via.d(i, j, k) * via.w(i, j, k) - dt_ / dy
                         * (flux_.w(i, j + 1, k) - flux_.w(i, j, k)))
                         / via_new->d(i, j, k);
        via_new->p(i, j, k) = via.p(i, j, k) + (gamma_ - 1.0) * (0.5 * (
                        via.d(i, j, k) * (via.u(i, j, k) * via.u(i, j, k)
                                      + via.v(i, j, k) * via.v(i, j, k)
                                      + via.w(i, j, k) * via.w(i, j, k))
                    - via_new->d(i, j, k) * (via_new->u(i, j, k) * via_new->u(i, j, k)
                                        + via_new->v(i, j, k) * via_new->v(i, j, k)
                                        + via_new->w(i, j, k) * via_new->w(i, j, k))))
                      - dt_ / dy * (flux_.p(i, j + 1, k) - flux_.p(i, j, k));
      }
    }
  }

  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        flux_.d(i, j, k) = 0.0;
        flux_.u(i, j, k) = 0.0;
        flux_.v(i, j, k) = 0.0;
        flux_.w(i, j, k) = 0.0;
        flux_.p(i, j, k) = 0.0;
      }
    }
  }
}

template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
            ::AdvCipRK3Y(Int32 l, Int32 i, Int32 j, Int32 k,
                         FluidQuantity3D rk[4],
                         const FluidQuantity3DArray &sia,
                         const FluidQuantity3DArray &via
                         ) noexcept
{
  Real lambda[3];
  Real dk[3], uk[3], vk[3], wk[3], pk[3];
  Real a[3] = {};

  switch (l) {
    default:
    case 0:
      a[0] = 1.0;
      a[1] = 0.0;
      a[2] = 0.0;
      break;
    case 1:
      a[0] = 0.25;
      a[1] = 0.25;
      a[2] = 0.0;
      break;
    case 2:
      a[0] = 1.0 / 6.0;
      a[1] = 1.0 / 6.0;
      a[2] = 2.0 / 3.0;
      break;
  }

  Real rM = a[0] * rk[0].d() + a[1] * rk[1].d() + a[2] * rk[2].d();
  Real vM = a[0] * rk[0].v() + a[1] * rk[1].v() + a[2] * rk[2].v();
  Real pM = a[0] * rk[0].p() + a[1] * rk[1].p() + a[2] * rk[2].p();
  Real cM = std::sqrt(gamma_ * pM / rM);
  assert(!std::isnan(cM));

  lambda[0] = -(vM     ) * dt_;
  lambda[1] = -(vM + cM) * dt_;
  lambda[2] = -(vM - cM) * dt_;

  //lambda 0
  dk[0] = cipy_(i, j, k, lambda[0], y_,
                   FQAccessor<TagD>(sia), FQAccessor<TagD>(via));
  uk[0] = cipy_(i, j, k, lambda[0], y_,
                   FQAccessor<TagU>(sia), FQAccessor<TagU>(via));
  wk[0] = cipy_(i, j, k, lambda[0], y_,
                   FQAccessor<TagW>(sia), FQAccessor<TagW>(via));
  pk[0] = cipy_(i, j, k, lambda[0], y_,
                   FQAccessor<TagP>(sia), FQAccessor<TagP>(via));

  //lambda 1
  vk[1] = cipy_(i, j, k, lambda[1], y_,
                   FQAccessor<TagV>(sia), FQAccessor<TagV>(via));
  pk[1] = cipy_(i, j, k, lambda[1], y_,
                   FQAccessor<TagP>(sia), FQAccessor<TagP>(via));

  //lambda 2
  vk[2] = cipy_(i, j, k, lambda[2], y_,
                   FQAccessor<TagV>(sia), FQAccessor<TagV>(via));
  pk[2] = cipy_(i, j, k, lambda[2], y_,
                   FQAccessor<TagP>(sia), FQAccessor<TagP>(via));

  //Characteristic_relationship
  rk[l + 1].u() = uk[0];
  rk[l + 1].v() = 0.5 * (vk[1] + vk[2] + (pk[1] - pk[2]) / (rM * cM));
  rk[l + 1].w() = wk[0];
  rk[l + 1].p() = 0.5 * (pk[1] + pk[2] + (vk[1] - vk[2]) * rM * cM);
  rk[l + 1].d() = dk[0] + (rk[l + 1].p() - pk[0]) / (cM * cM);
}
