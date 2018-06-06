#ifndef FLUID1D_H
#define FLUID1D_H

#include <functional>
#include <string>

#include "global_definitions.h"

class Fluid1D
{
public:
  Fluid1D() = default;
  Fluid1D(Int32 size, Real xmin, Real xmax, Real gamma,
          std::function<void(const Rvec&,
                             Rvec*,
                             Rvec*,
                             Rvec*)> init_func);

  void MainLoop(Int32 loop_max, Int32 out_step, const std::string &fname_base);
  void Update() noexcept;
  void Output(const std::string &filename, Int32 step) const;

  Fluid1D& Initialize(Int32 size, Real xmin, Real xmax, Real gamma,
                      std::function<void(const Rvec&,
                                         Rvec*,
                                         Rvec*,
                                         Rvec*)> init_func);
private:
  Int32 size_;
  Real xmin_, xmax_;
  Rvec x_;
  Real dt_;
  Real time_;

  Real gamma_{1.4};

  Rvec density_prim_sia_;
  Rvec velocity_prim_sia_;
  Rvec pressure_prim_sia_;

  Rvec density_prim_sia_new_;
  Rvec velocity_prim_sia_new_;
  Rvec pressure_prim_sia_new_;

  Rvec density_prim_via_;
  Rvec velocity_prim_via_;
  Rvec pressure_prim_via_;

  Rvec density_prim_via_new_;
  Rvec velocity_prim_via_new_;
  Rvec pressure_prim_via_new_;

  Rvec density_cons_via_;
  Rvec momentum_cons_via_;
  Rvec energy_cons_via_;

  Rvec density_flux_;
  Rvec momentum_flux_;
  Rvec energy_flux_;

  void NewDt() noexcept;

  void BoundaryConditions() noexcept;

  void Advection() noexcept;
  void AdvectionSia() noexcept;
  void AdvectionVia() noexcept;

  void AdvCipRK3(Int32 l, Int32 i,
                 Real dens[4], Real velo[4], Real pres[4]) noexcept;

  Real CipCsl3(Int32 i, Real v, const Rvec& sia, const Rvec& via) noexcept;
  Real CipCsl3(Int32 i, Real v, Real v_bck, Real v_fwd,
               const Rvec& sia, const Rvec& via) noexcept;

  void UpdateValue() noexcept;
};

#endif // FLUID1D_H
