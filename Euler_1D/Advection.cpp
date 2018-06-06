#include <cmath>

#include "global_definitions.h"
#include "Fluid1D.h"


void Fluid1D::Advection() noexcept
{
  AdvectionSia();
  AdvectionVia();
}

/*
void Fluid1D::AdvectionSia() noexcept
{
  auto &dps = density_prim_sia_;
  auto &vps = velocity_prim_sia_;
  auto &pps = pressure_prim_sia_;
  auto &dpsn = density_prim_sia_new_;
  auto &vpsn = velocity_prim_sia_new_;
  auto &ppsn = pressure_prim_sia_new_;

  Real dens[3], velo[3], pres[3];

  for(auto i = GHOST; i < size_ + 1 - GHOST; ++i){
    auto cs = std::sqrt(gamma_ * pps[i] / dps[i]);
    auto cs_bck = std::sqrt(gamma_ * pps[i - 1] / dps[i - 1]);
    auto cs_fwd = std::sqrt(gamma_ * pps[i + 1] / dps[i + 1]);
    auto v = vps[i];
    auto v_bck = vps[i - 1];
    auto v_fwd = vps[i + 1];


    dens[0] = CipCsl3(i, v, v_bck, v_fwd, dps, density_prim_via_);
    pres[0] = CipCsl3(i, v, v_bck, v_fwd, pps, pressure_prim_via_);

    velo[1] = CipCsl3(i, v + cs, v_bck + cs_bck, v_fwd + cs_fwd,
                      vps, velocity_prim_via_);
    pres[1] = CipCsl3(i, v + cs, v_bck + cs_bck, v_fwd + cs_fwd,
                      pps, pressure_prim_via_);

    velo[2] = CipCsl3(i, v - cs, v_bck - cs_bck, v_fwd - cs_fwd,
                      vps, velocity_prim_via_);
    pres[2] = CipCsl3(i, v - cs, v_bck - cs_bck, v_fwd - cs_fwd,
                      pps, pressure_prim_via_);

    dens[0] = CipCsl3(i, v, dps, density_prim_via_);
    pres[0] = CipCsl3(i, v, pps, pressure_prim_via_);

    velo[1] = CipCsl3(i, v + cs, vps, velocity_prim_via_);
    pres[1] = CipCsl3(i, v + cs, pps, pressure_prim_via_);

    velo[2] = CipCsl3(i, v - cs, vps, velocity_prim_via_);
    pres[2] = CipCsl3(i, v - cs, pps, pressure_prim_via_);

    ppsn[i] = 0.5 * (pres[1] + pres[2] + dps[i] * cs * (velo[1] - velo[2]));
    vpsn[i] = 0.5 * (velo[1] + velo[2] + (pres[1] - pres[2]) / (dps[i] * cs));
    dpsn[i] = dens[0] + (ppsn[i] - pres[0]) / (cs * cs);

    density_flux_[i] = dpsn[i] * vpsn[i];
    momentum_flux_[i] = dpsn[i] * vpsn[i] * vpsn[i] + ppsn[i];
    auto energy = ppsn[i] / (gamma_ - 1.0) + 0.5 * dpsn[i] * vpsn[i] * vpsn[i];
    energy_flux_[i] = vpsn[i] * (energy + ppsn[i]);
  }
}
*/

void Fluid1D::AdvectionSia() noexcept
{
  Real dens[4], velo[4], pres[4], engy[3];
  Real dflux[3], mflux[3], eflux[3];

  for(auto i = GHOST; i < size_ + 1 - GHOST; ++i){
    dens[0] = density_prim_sia_[i];
    velo[0] = velocity_prim_sia_[i];
    pres[0] = pressure_prim_sia_[i];

    AdvCipRK3(0, i, dens, velo, pres);
    AdvCipRK3(1, i, dens, velo, pres);
    AdvCipRK3(2, i, dens, velo, pres);

    for(auto l = 0; l < 3; ++l){
      engy[l] = pres[l] / (gamma_ - 1.0) + 0.5 * dens[l] * velo[l] * velo[l];
    }

    for(auto l = 0; l < 3; ++l){
      dflux[l] = dens[l] * velo[l];
      mflux[l] = dens[l] * velo[l] * velo[l] + pres[l];
      eflux[l] = velo[l] * (engy[l] + pres[l]);
    }

    density_prim_sia_new_[i] = dens[3];
    velocity_prim_sia_new_[i] = velo[3];
    pressure_prim_sia_new_[i] = pres[3];

    density_flux_[i] = (dflux[0] + dflux[1] + 4.0 * dflux[2]) / 6.0;
    momentum_flux_[i] = (mflux[0] + mflux[1] + 4.0 * mflux[2]) / 6.0;
    energy_flux_[i] = (eflux[0] + eflux[1] + 4.0 * eflux[2]) / 6.0;
  }
}

void Fluid1D::AdvCipRK3(Int32 l, Int32 i,
                        Real dens[4], Real velo[4], Real pres[4]) noexcept
{
  Real lambda[3];
  Real rk[3],vk[3],pk[3];
  Real a[3] = {0.0, 0.0, 0.0};

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
  Real vM = a[0] * velo[0] + a[1] * velo[1] + a[2] * velo[2];
  Real pM = a[0] * pres[0] + a[1] * pres[1] + a[2] * pres[2];
  Real cM = std::sqrt(gamma_ * pM / rM);

  lambda[0] = vM;
  lambda[1] = vM + cM;
  lambda[2] = vM - cM;

  //lambda 0
  rk[0] = CipCsl3(i, lambda[0], density_prim_sia_, density_prim_via_);
  pk[0] = CipCsl3(i, lambda[0], pressure_prim_sia_, pressure_prim_via_);

  //lambda 1
  vk[1] = CipCsl3(i, lambda[1], velocity_prim_sia_, velocity_prim_via_);
  pk[1] = CipCsl3(i, lambda[1], pressure_prim_sia_, pressure_prim_via_);

  //lambda 2
  vk[2] = CipCsl3(i, lambda[2], velocity_prim_sia_, velocity_prim_via_);
  pk[2] = CipCsl3(i, lambda[2], pressure_prim_sia_, pressure_prim_via_);

  //Characteristic_relationship
  velo[l + 1] = 0.5 * (vk[1] + vk[2] + (pk[1] - pk[2]) / (rM * cM));
  pres[l + 1] = 0.5 * (pk[1] + pk[2] + (vk[1] - vk[2]) * rM * cM);
  dens[l + 1] = rk[0] + (pres[l + 1] - pk[0]) / (cM * cM);
}

void Fluid1D::AdvectionVia() noexcept
{
  for(auto i = GHOST; i < size_ - GHOST; ++i){
    auto dx = x_[i + 1] - x_[i];
    density_cons_via_[i] -= dt_ / dx
                          * (density_flux_[i + 1] - density_flux_[i]);
    momentum_cons_via_[i] -= dt_ / dx
                           * (momentum_flux_[i + 1] - momentum_flux_[i]);
    energy_cons_via_[i] -= dt_ / dx
                         * (energy_flux_[i + 1] - energy_flux_[i]);

    density_prim_via_new_[i] = density_cons_via_[i];
    velocity_prim_via_new_[i] = momentum_cons_via_[i] / density_cons_via_[i];
    pressure_prim_via_new_[i] = (gamma_ - 1.0) * (energy_cons_via_[i]
                              - 0.5 * momentum_cons_via_[i]
                              * momentum_cons_via_[i] / density_cons_via_[i]);
  }
}
