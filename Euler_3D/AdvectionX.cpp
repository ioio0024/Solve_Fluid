#include <cassert>
#include <cmath>

#include "global_definitions.h"
#include "Fluid3D.h"

void Fluid3D::AdvectionX(Int32 ny, Int32 nz,
                         const Vector3DR &d_p,
                         const Vector3DR &u_p,
                         const Vector3DR &v_p,
                         const Vector3DR &w_p,
                         const Vector3DR &p_p,
                         const Vector3DR &d_v,
                         const Vector3DR &u_v,
                         const Vector3DR &v_v,
                         const Vector3DR &w_v,
                         const Vector3DR &p_v,
                         Vector3DR *d_pn,
                         Vector3DR *u_pn,
                         Vector3DR *v_pn,
                         Vector3DR *w_pn,
                         Vector3DR *p_pn,
                         Vector3DR *d_vn,
                         Vector3DR *u_vn,
                         Vector3DR *v_vn,
                         Vector3DR *w_vn,
                         Vector3DR *p_vn) noexcept
{
  Real dens[4], velx[4], vely[4], velz[4], pres[4], engy[3];
  Real dflux[3], lflux[3], mflux[3], nflux[3], eflux[3];

  for(auto k = GHOST; k < nz - GHOST; ++k){
    for(auto j = GHOST; j < ny - GHOST; ++j){
      for(auto i = GHOST; i < nx_ + 1 - GHOST; ++i){
        dens[0] = d_p(i, j, k);
        velx[0] = u_p(i, j, k);
        vely[0] = v_p(i, j, k);
        velz[0] = w_p(i, j, k);
        pres[0] = p_p(i, j, k);

        AdvCipRK3X(0, i, j, k, dens, velx, vely, velz, pres,
                   d_p, u_p, v_p, w_p, p_p,
                   d_v, u_v, v_v, w_v, p_v);
        AdvCipRK3X(1, i, j, k, dens, velx, vely, velz, pres,
                   d_p, u_p, v_p, w_p, p_p,
                   d_v, u_v, v_v, w_v, p_v);
        AdvCipRK3X(2, i, j, k, dens, velx, vely, velz, pres,
                   d_p, u_p, v_p, w_p, p_p,
                   d_v, u_v, v_v, w_v, p_v);

        for(auto l = 0; l < 3; ++l){
          engy[l] = pres[l] + 0.5 * (gamma_ - 1.0)
                  * dens[l] * (velx[l] * velx[l]
                             + vely[l] * vely[l]
                             + velz[l] * velz[l]);
        }

        for(auto l = 0; l < 3; ++l){
          dflux[l] = dens[l] * velx[l];
          lflux[l] = dens[l] * velx[l] * velx[l] + pres[l];
          mflux[l] = dens[l] * velx[l] * vely[l];
          nflux[l] = dens[l] * velx[l] * velz[l];
          eflux[l] = velx[l] * (engy[l] + (gamma_ - 1.0) * pres[l]);
        }

        (*d_pn)(i, j, k) = dens[3];
        (*u_pn)(i, j, k) = velx[3];
        (*v_pn)(i, j, k) = vely[3];
        (*w_pn)(i, j, k) = velz[3];
        (*p_pn)(i, j, k) = pres[3];

        dflux_(i, j, k) = (dflux[0] + dflux[1] + 4.0 * dflux[2]) / 6.0;
        lflux_(i, j, k) = (lflux[0] + lflux[1] + 4.0 * lflux[2]) / 6.0;
        mflux_(i, j, k) = (mflux[0] + mflux[1] + 4.0 * mflux[2]) / 6.0;
        nflux_(i, j, k) = (nflux[0] + nflux[1] + 4.0 * nflux[2]) / 6.0;
        eflux_(i, j, k) = (eflux[0] + eflux[1] + 4.0 * eflux[2]) / 6.0;
      }
    }
  }

  // VIA
  for(auto k = GHOST; k < nz - GHOST; ++k){
    for(auto j = GHOST; j < ny - GHOST; ++j){
      for(auto i = GHOST; i < nx_ - GHOST; ++i){
        auto dx = x_[i + 1] - x_[i];
        (*d_vn)(i, j, k) = d_v(i, j, k) - dt_ / dx
            * (dflux_(i + 1, j, k) - dflux_(i, j, k));
        (*u_vn)(i, j, k) = (d_v(i, j, k) * u_v(i, j, k) - dt_ / dx
                         * (lflux_(i + 1, j, k) - lflux_(i, j, k)))
                         / (*d_vn)(i, j, k);
        (*v_vn)(i, j, k) = (d_v(i, j, k) * v_v(i, j, k) - dt_ / dx
                         * (mflux_(i + 1, j, k) - mflux_(i, j, k)))
                         / (*d_vn)(i, j, k);
        (*w_vn)(i, j, k) = (d_v(i, j, k) * w_v(i, j, k) - dt_ / dx
                         * (nflux_(i + 1, j, k) - nflux_(i, j, k)))
                         / (*d_vn)(i, j, k);
        (*p_vn)(i, j, k) = p_v(i, j, k) + (gamma_ - 1.0) * (0.5 * (
                        d_v(i, j, k) * (u_v(i, j, k) * u_v(i, j, k)
                                      + v_v(i, j, k) * v_v(i, j, k)
                                      + w_v(i, j, k) * w_v(i, j, k))
                    - (*d_vn)(i, j, k) * ((*u_vn)(i, j, k) * (*u_vn)(i, j, k)
                                        + (*v_vn)(i, j, k) * (*v_vn)(i, j, k)
                                        + (*w_vn)(i, j, k) * (*w_vn)(i, j, k))))
                      - dt_ / dx * (eflux_(i + 1, j, k) - eflux_(i, j, k));
      }
    }
  }

  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        dflux_(i, j, k) = 0.0;
        lflux_(i, j, k) = 0.0;
        mflux_(i, j, k) = 0.0;
        nflux_(i, j, k) = 0.0;
        eflux_(i, j, k) = 0.0;
      }
    }
  }
}

void Fluid3D::AdvCipRK3X(Int32 l, Int32 i, Int32 j, Int32 k,
                         Real dens[4], Real velx[4], Real vely[4],
                         Real velz[4], Real pres[4],
                         const Vector3DR &d_p, const Vector3DR &u_p,
                         const Vector3DR &v_p, const Vector3DR &w_p,
                         const Vector3DR &p_p,
                         const Vector3DR &d_v, const Vector3DR &u_v,
                         const Vector3DR &v_v, const Vector3DR &w_v,
                         const Vector3DR &p_v) noexcept
{
  Real lambda[3];
  Real rk[3], uk[3], vk[3], wk[3], pk[3];
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
  rk[0] = CipCsl3X(i, j, k, lambda[0], d_p, d_v);
  vk[0] = CipCsl3X(i, j, k, lambda[0], v_p, v_v);
  wk[0] = CipCsl3X(i, j, k, lambda[0], w_p, w_v);
  pk[0] = CipCsl3X(i, j, k, lambda[0], p_p, p_v);

  //lambda 1
  uk[1] = CipCsl3X(i, j, k, lambda[1], u_p, u_v);
  pk[1] = CipCsl3X(i, j, k, lambda[1], p_p, p_v);

  //lambda 2
  uk[2] = CipCsl3X(i, j, k, lambda[2], u_p, u_v);
  pk[2] = CipCsl3X(i, j, k, lambda[2], p_p, p_v);

  //Characteristic_relationship
  velx[l + 1] = 0.5 * (uk[1] + uk[2] + (pk[1] - pk[2]) / (rM * cM));
  vely[l + 1] = vk[0];
  velz[l + 1] = wk[0];
  pres[l + 1] = 0.5 * (pk[1] + pk[2] + (uk[1] - uk[2]) * rM * cM);
  dens[l + 1] = rk[0] + (pres[l + 1] - pk[0]) / (cM * cM);
}
