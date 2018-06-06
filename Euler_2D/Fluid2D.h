#ifndef FLUID2D_H
#define FLUID2D_H

#include <functional>
#include <string>

#include "global_definitions.h"
#include "Vector2D.h"

class Fluid2D
{
public:
  Fluid2D() = default;
  Fluid2D(Int32 size_x, Int32 size_y,
          Real xmin, Real xmax,
          Real ymin, Real ymax,
          Real gamma,
          std::function<void(const Rvec&, const Rvec&,
                             Vector2D*,
                             Vector2D*,
                             Vector2D*,
                             Vector2D*)> init_func)
  {
    Initialize(size_x, size_y, xmin, xmax, ymin, ymax, gamma, init_func);
  }

  void MainLoop(Int32 loop_max, Int32 out_step, const std::string &fname_base);
  void Update() noexcept;
  void Output(const std::string &filename, Int32 step) const;
  void OutputSiaX(const std::string &filename, Int32 step) const;
  void OutputSiaY(const std::string &filename, Int32 step) const;

  Fluid2D& Initialize(Int32 size_x, Int32 size_y,
                      Real xmin, Real xmax,
                      Real ymin, Real ymax,
                      Real gamma,
                      std::function<void(const Rvec&, const Rvec&,
                                         Vector2D*,
                                         Vector2D*,
                                         Vector2D*,
                                         Vector2D*)> init_func);
private:
  Int32 size_x_, size_y_;
  Real xmin_, xmax_;
  Real ymin_, ymax_;
  Rvec x_;
  Rvec y_;
  Real dt_;
  Real time_;

  Real gamma_{1.4};

  Vector2D density_prim_pv_;
  Vector2D vx_prim_pv_;
  Vector2D vy_prim_pv_;
  Vector2D pressure_prim_pv_;

  Vector2D density_prim_sia_x_;
  Vector2D vx_prim_sia_x_;
  Vector2D vy_prim_sia_x_;
  Vector2D pressure_prim_sia_x_;

  Vector2D density_prim_sia_y_;
  Vector2D vx_prim_sia_y_;
  Vector2D vy_prim_sia_y_;
  Vector2D pressure_prim_sia_y_;

  Vector2D density_prim_via_;
  Vector2D vx_prim_via_;
  Vector2D vy_prim_via_;
  Vector2D pressure_prim_via_;

  Vector2D density_prim_pv_new_;
  Vector2D vx_prim_pv_new_;
  Vector2D vy_prim_pv_new_;
  Vector2D pressure_prim_pv_new_;

  Vector2D density_prim_sia_x_new_;
  Vector2D vx_prim_sia_x_new_;
  Vector2D vy_prim_sia_x_new_;
  Vector2D pressure_prim_sia_x_new_;

  Vector2D density_prim_sia_y_new_;
  Vector2D vx_prim_sia_y_new_;
  Vector2D vy_prim_sia_y_new_;
  Vector2D pressure_prim_sia_y_new_;

  Vector2D density_prim_via_new_;
  Vector2D vx_prim_via_new_;
  Vector2D vy_prim_via_new_;
  Vector2D pressure_prim_via_new_;

  Vector2D dflux_x_;
  Vector2D lflux_x_;
  Vector2D mflux_x_;
  Vector2D eflux_x_;

  Vector2D dflux_y_;
  Vector2D lflux_y_;
  Vector2D mflux_y_;
  Vector2D eflux_y_;

  void NewDt() noexcept;

  void BoundaryConditions() noexcept;
  void BoundaryConditionsX() noexcept;
  void BoundaryConditionsY() noexcept;

  void Advection() noexcept;
  void AdvectionX(Int32 n,
                  const Vector2D& d_p, const Vector2D& u_p,
                  const Vector2D& v_p, const Vector2D& p_p,
                  const Vector2D& d_v, const Vector2D& u_v,
                  const Vector2D& v_v, const Vector2D& p_v,
                  Vector2D* d_pn, Vector2D* u_pn,
                  Vector2D* v_pn, Vector2D* p_pn,
                  Vector2D* d_vn, Vector2D* u_vn,
                  Vector2D* v_vn,Vector2D* p_vn) noexcept;
  void AdvectionY(Int32 n,
                  const Vector2D& d_p, const Vector2D& u_p,
                  const Vector2D& v_p, const Vector2D& p_p,
                  const Vector2D& d_v, const Vector2D& u_v,
                  const Vector2D& v_v, const Vector2D& p_v,
                  Vector2D* d_pn, Vector2D* u_pn,
                  Vector2D* v_pn, Vector2D* p_pn,
                  Vector2D* d_vn, Vector2D* u_vn,
                  Vector2D* v_vn, Vector2D* p_vn) noexcept;

  void AdvCipRK3X(Int32 l, Int32 i, Int32 j,
                  Real dens[4], Real velx[4], Real vely[4], Real pres[4],
                  const Vector2D& d_p, const Vector2D& u_p,
                  const Vector2D& v_p, const Vector2D& p_p,
                  const Vector2D& d_v, const Vector2D& u_v,
                  const Vector2D& v_v, const Vector2D& p_v) noexcept;

  void AdvCipRK3Y(Int32 l, Int32 i, Int32 j,
                  Real dens[4], Real velx[4], Real vely[4], Real pres[4],
                  const Vector2D& d_p, const Vector2D& u_p,
                  const Vector2D& v_p, const Vector2D& p_p,
                  const Vector2D& d_v, const Vector2D& u_v,
                  const Vector2D& v_v, const Vector2D& p_v) noexcept;

  Real CipCsl3X(Int32 i, Int32 j, Real v,
                const Vector2D& sia, const Vector2D& via) noexcept;
  Real CipCsl3Y(Int32 i, Int32 j, Real v,
                const Vector2D& sia, const Vector2D& via) noexcept;

  void UpdateValue() noexcept;
};

#endif // FLUID2D_H
