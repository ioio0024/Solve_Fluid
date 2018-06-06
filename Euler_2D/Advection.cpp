#include <cassert>
#include <cmath>

#include "global_definitions.h"
#include "Fluid2D.h"

void Fluid2D::Advection() noexcept
{
  AdvectionX(size_y_ + 1,
             density_prim_pv_, vx_prim_pv_,
             vy_prim_pv_, pressure_prim_pv_,
             density_prim_sia_y_, vx_prim_sia_y_,
             vy_prim_sia_y_, pressure_prim_sia_y_,
             &density_prim_pv_new_, &vx_prim_pv_new_,
             &vy_prim_pv_new_, &pressure_prim_pv_new_,
             &density_prim_sia_y_new_, &vx_prim_sia_y_new_,
             &vy_prim_sia_y_new_, &pressure_prim_sia_y_new_);
  AdvectionX(size_y_,
             density_prim_sia_x_, vx_prim_sia_x_,
             vy_prim_sia_x_, pressure_prim_sia_x_,
             density_prim_via_, vx_prim_via_,
             vy_prim_via_, pressure_prim_via_,
             &density_prim_sia_x_new_, &vx_prim_sia_x_new_,
             &vy_prim_sia_x_new_, &pressure_prim_sia_x_new_,
             &density_prim_via_new_, &vx_prim_via_new_,
             &vy_prim_via_new_, &pressure_prim_via_new_);
  UpdateValue();
  BoundaryConditions();

  AdvectionY(size_x_ + 1,
             density_prim_pv_, vx_prim_pv_,
             vy_prim_pv_, pressure_prim_pv_,
             density_prim_sia_x_, vx_prim_sia_x_,
             vy_prim_sia_x_, pressure_prim_sia_x_,
             &density_prim_pv_new_, &vx_prim_pv_new_,
             &vy_prim_pv_new_, &pressure_prim_pv_new_,
             &density_prim_sia_x_new_, &vx_prim_sia_x_new_,
             &vy_prim_sia_x_new_, &pressure_prim_sia_x_new_);
  AdvectionY(size_x_,
             density_prim_sia_y_, vx_prim_sia_y_,
             vy_prim_sia_y_, pressure_prim_sia_y_,
             density_prim_via_, vx_prim_via_,
             vy_prim_via_, pressure_prim_via_,
             &density_prim_sia_y_new_, &vx_prim_sia_y_new_,
             &vy_prim_sia_y_new_, &pressure_prim_sia_y_new_,
             &density_prim_via_new_, &vx_prim_via_new_,
             &vy_prim_via_new_, &pressure_prim_via_new_);
  UpdateValue();
}

void Fluid2D::AdvectionX(Int32 n,
                         const Vector2D &d_p, const Vector2D &u_p,
                         const Vector2D &v_p, const Vector2D &p_p,
                         const Vector2D &d_v, const Vector2D &u_v,
                         const Vector2D &v_v, const Vector2D &p_v,
                         Vector2D *d_pn, Vector2D *u_pn,
                         Vector2D *v_pn, Vector2D *p_pn,
                         Vector2D *d_vn, Vector2D *u_vn,
                         Vector2D *v_vn, Vector2D *p_vn) noexcept
{
  Real dens[4], velx[4], vely[4], pres[4], engy[3];
  Real dflux[3], lflux[3], mflux[3], eflux[3];

  for(auto j = GHOST; j < n - GHOST; ++j){
    for(auto i = GHOST; i < size_x_ + 1 - GHOST; ++i){
      dens[0] = d_p(i, j);
      velx[0] = u_p(i, j);
      vely[0] = v_p(i, j);
      pres[0] = p_p(i, j);

      AdvCipRK3X(0, i, j, dens, velx, vely, pres,
                 d_p, u_p, v_p, p_p,
                 d_v, u_v, v_v, p_v);
      AdvCipRK3X(1, i, j, dens, velx, vely, pres,
                 d_p, u_p, v_p, p_p,
                 d_v, u_v, v_v, p_v);
      AdvCipRK3X(2, i, j, dens, velx, vely, pres,
                 d_p, u_p, v_p, p_p,
                 d_v, u_v, v_v, p_v);

      for(auto l = 0; l < 3; ++l){
        engy[l] = pres[l] + 0.5 * (gamma_ - 1.0)
                * dens[l] * (velx[l] * velx[l] + vely[l] * vely[l]);
      }

      for(auto l = 0; l < 3; ++l){
        dflux[l] = dens[l] * velx[l];
        lflux[l] = dens[l] * velx[l] * velx[l] + pres[l];
        mflux[l] = dens[l] * velx[l] * vely[l];
        eflux[l] = velx[l] * (engy[l] + (gamma_ - 1.0) * pres[l]);
      }

      (*d_pn)(i, j) = dens[3];
      (*u_pn)(i, j) = velx[3];
      (*v_pn)(i, j) = vely[3];
      (*p_pn)(i, j) = pres[3];

      dflux_x_(i, j) = (dflux[0] + dflux[1] + 4.0 * dflux[2]) / 6.0;
      lflux_x_(i, j) = (lflux[0] + lflux[1] + 4.0 * lflux[2]) / 6.0;
      mflux_x_(i, j) = (mflux[0] + mflux[1] + 4.0 * mflux[2]) / 6.0;
      eflux_x_(i, j) = (eflux[0] + eflux[1] + 4.0 * eflux[2]) / 6.0;
    }
  }

  for(auto j = GHOST; j < n - GHOST; ++j){
    for(auto i = GHOST; i < size_x_ - GHOST; ++i){
      auto dx = x_[i + 1] - x_[i];
      (*d_vn)(i, j) = d_v(i, j) - dt_ / dx
                    * (dflux_x_(i + 1, j) - dflux_x_(i, j));
      (*u_vn)(i, j) = (d_v(i, j) * u_v(i, j) - dt_ / dx
                    * (lflux_x_(i + 1, j) - lflux_x_(i, j)))
                    / (*d_vn)(i, j);
      (*v_vn)(i, j) = (d_v(i, j) * v_v(i, j) - dt_ / dx
                    * (mflux_x_(i + 1, j) - mflux_x_(i, j)))
                    / (*d_vn)(i, j);
      (*p_vn)(i, j) = p_v(i, j) + (gamma_ - 1.0) * (0.5 * (
                      d_v(i, j) * (u_v(i, j) * u_v(i, j) + v_v(i, j) * v_v(i, j))
                    - (*d_vn)(i, j) * ((*u_vn)(i, j) * (*u_vn)(i, j)
                                     + (*v_vn)(i, j) * (*v_vn)(i, j))))
                    - dt_ / dx * (eflux_x_(i + 1, j) - eflux_x_(i, j));
    }
  }

  for(auto j = 0; j < size_y_ + 1; ++j){
    for(auto i = 0; i < size_x_ + 1; ++i){
      dflux_x_(i, j) = 0.0;
      lflux_x_(i, j) = 0.0;
      mflux_x_(i, j) = 0.0;
      eflux_x_(i, j) = 0.0;
    }
  }
}

void Fluid2D::AdvectionY(Int32 n,
                         const Vector2D &d_p, const Vector2D &u_p,
                         const Vector2D &v_p, const Vector2D &p_p,
                         const Vector2D &d_v, const Vector2D &u_v,
                         const Vector2D &v_v, const Vector2D &p_v,
                         Vector2D *d_pn, Vector2D *u_pn,
                         Vector2D *v_pn, Vector2D *p_pn,
                         Vector2D *d_vn, Vector2D *u_vn,
                         Vector2D *v_vn, Vector2D *p_vn) noexcept
{
  Real dens[4], velx[4], vely[4], pres[4], engy[3];
  Real dflux[3], lflux[3], mflux[3], eflux[3];

  for(auto j = GHOST; j < size_y_ + 1 - GHOST; ++j){
    for(auto i = GHOST; i < n - GHOST; ++i){
      dens[0] = d_p(i, j);
      velx[0] = u_p(i, j);
      vely[0] = v_p(i, j);
      pres[0] = p_p(i, j);

      AdvCipRK3Y(0, i, j, dens, velx, vely, pres,
                 d_p, u_p, v_p, p_p,
                 d_v, u_v, v_v, p_v);
      AdvCipRK3Y(1, i, j, dens, velx, vely, pres,
                 d_p, u_p, v_p, p_p,
                 d_v, u_v, v_v, p_v);
      AdvCipRK3Y(2, i, j, dens, velx, vely, pres,
                 d_p, u_p, v_p, p_p,
                 d_v, u_v, v_v, p_v);

      for(auto l = 0; l < 3; ++l){
        engy[l] = pres[l] + 0.5 * (gamma_ - 1.0)
                * dens[l] * (velx[l] * velx[l] + vely[l] * vely[l]);
      }

      for(auto l = 0; l < 3; ++l){
        dflux[l] = dens[l] * vely[l];
        lflux[l] = dens[l] * velx[l] * vely[l];
        mflux[l] = dens[l] * vely[l] * vely[l] + pres[l];
        eflux[l] = vely[l] * (engy[l] + (gamma_ - 1.0) * pres[l]);
      }

      (*d_pn)(i, j) = dens[3];
      (*u_pn)(i, j) = velx[3];
      (*v_pn)(i, j) = vely[3];
      (*p_pn)(i, j) = pres[3];

      dflux_y_(i, j) = (dflux[0] + dflux[1] + 4.0 * dflux[2]) / 6.0;
      lflux_y_(i, j) = (lflux[0] + lflux[1] + 4.0 * lflux[2]) / 6.0;
      mflux_y_(i, j) = (mflux[0] + mflux[1] + 4.0 * mflux[2]) / 6.0;
      eflux_y_(i, j) = (eflux[0] + eflux[1] + 4.0 * eflux[2]) / 6.0;
    }
  }

  for(auto j = GHOST; j < size_y_ - GHOST; ++j){
    for(auto i = GHOST; i < n - GHOST; ++i){
      auto dy = y_[j + 1] - y_[j];
      (*d_vn)(i, j) = d_v(i, j) - dt_ / dy
                    * (dflux_y_(i, j + 1) - dflux_y_(i, j));
      (*u_vn)(i, j) = (d_v(i, j) * u_v(i, j) - dt_ / dy
                    * (lflux_y_(i, j + 1) - lflux_y_(i, j)))
                    / (*d_vn)(i, j);
      (*v_vn)(i, j) = (d_v(i, j) * v_v(i, j) - dt_ / dy
                    * (mflux_y_(i, j + 1) - mflux_y_(i, j)))
                    / (*d_vn)(i, j);
      (*p_vn)(i, j) = p_v(i, j) + (gamma_ - 1.0) * (0.5 * (
                      d_v(i, j) * (u_v(i, j) * u_v(i, j) + v_v(i, j) * v_v(i, j))
                    - (*d_vn)(i, j) * ((*u_vn)(i, j) * (*u_vn)(i, j)
                                     + (*v_vn)(i, j) * (*v_vn)(i, j))))
                    - dt_ / dy * (eflux_y_(i, j + 1) - eflux_y_(i, j));
    }
  }

  for(auto j = 0; j < size_y_ + 1; ++j){
    for(auto i = 0; i < size_x_ + 1; ++i){
      dflux_y_(i, j) = 0.0;
      lflux_y_(i, j) = 0.0;
      mflux_y_(i, j) = 0.0;
      eflux_y_(i, j) = 0.0;
    }
  }
}

void Fluid2D::AdvCipRK3X(Int32 l, Int32 i, Int32 j,
                        Real dens[4], Real velx[4], Real vely[4],Real pres[4],
                        const Vector2D& d_p, const Vector2D& u_p,
                        const Vector2D& v_p, const Vector2D& p_p,
                        const Vector2D& d_v, const Vector2D& u_v,
                        const Vector2D& v_v, const Vector2D& p_v) noexcept
{
  Real lambda[3];
  Real rk[3],uk[3], vk[3],pk[3];
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

  Real rM = a[0] * dens[0] + a[1] * dens[1] + a[2] * dens[2];
  Real uM = a[0] * velx[0] + a[1] * velx[1] + a[2] * velx[2];
  Real pM = a[0] * pres[0] + a[1] * pres[1] + a[2] * pres[2];
  Real cM = std::sqrt(gamma_ * pM / rM);
  assert(!std::isnan(cM));

  lambda[0] = uM;
  lambda[1] = uM + cM;
  lambda[2] = uM - cM;

  //lambda 0
  rk[0] = CipCsl3X(i, j, lambda[0], d_p, d_v);
  vk[0] = CipCsl3X(i, j, lambda[0], v_p, v_v);
  pk[0] = CipCsl3X(i, j, lambda[0], p_p, p_v);

  //lambda 1
  uk[1] = CipCsl3X(i, j, lambda[1], u_p, u_v);
  pk[1] = CipCsl3X(i, j, lambda[1], p_p, p_v);

  //lambda 2
  uk[2] = CipCsl3X(i, j, lambda[2], u_p, u_v);
  pk[2] = CipCsl3X(i, j, lambda[2], p_p, p_v);

  //Characteristic_relationship
  velx[l + 1] = 0.5 * (uk[1] + uk[2] + (pk[1] - pk[2]) / (rM * cM));
  vely[l + 1] = vk[0];
  pres[l + 1] = 0.5 * (pk[1] + pk[2] + (uk[1] - uk[2]) * rM * cM);
  dens[l + 1] = rk[0] + (pres[l + 1] - pk[0]) / (cM * cM);
}

void Fluid2D::AdvCipRK3Y(Int32 l, Int32 i, Int32 j,
                        Real dens[4], Real velx[4], Real vely[4],Real pres[4],
                        const Vector2D& d_p, const Vector2D& u_p,
                        const Vector2D& v_p, const Vector2D& p_p,
                        const Vector2D& d_v, const Vector2D& u_v,
                        const Vector2D& v_v, const Vector2D& p_v) noexcept
{
  Real lambda[3];
  Real rk[3],uk[3], vk[3],pk[3];
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

  Real rM = a[0] * dens[0] + a[1] * dens[1] + a[2] * dens[2];
  Real vM = a[0] * vely[0] + a[1] * vely[1] + a[2] * vely[2];
  Real pM = a[0] * pres[0] + a[1] * pres[1] + a[2] * pres[2];
  Real cM = std::sqrt(gamma_ * pM / rM);
  assert(!std::isnan(cM));

  lambda[0] = vM;
  lambda[1] = vM + cM;
  lambda[2] = vM - cM;

  //lambda 0
  rk[0] = CipCsl3Y(i, j, lambda[0], d_p, d_v);
  uk[0] = CipCsl3Y(i, j, lambda[0], u_p, u_v);
  pk[0] = CipCsl3Y(i, j, lambda[0], p_p, p_v);

  //lambda 1
  vk[1] = CipCsl3Y(i, j, lambda[1], v_p, v_v);
  pk[1] = CipCsl3Y(i, j, lambda[1], p_p, p_v);

  //lambda 2
  vk[2] = CipCsl3Y(i, j, lambda[2], v_p, v_v);
  pk[2] = CipCsl3Y(i, j, lambda[2], p_p, p_v);

  //Characteristic_relationship
  velx[l + 1] = uk[0];
  vely[l + 1] = 0.5 * (vk[1] + vk[2] + (pk[1] - pk[2]) / (rM * cM));
  pres[l + 1] = 0.5 * (pk[1] + pk[2] + (vk[1] - vk[2]) * rM * cM);
  dens[l + 1] = rk[0] + (pres[l + 1] - pk[0]) / (cM * cM);
}
