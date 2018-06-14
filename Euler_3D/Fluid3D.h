#ifndef FLUID3D_H
#define FLUID3D_H

#include <functional>
#include <string>

#include "global_definitions.h"
#include "Vector3D.h"

class Fluid3D
{
public:
  Fluid3D() = default;
  Fluid3D(Int32 size_x, Int32 size_y, Int32 size_z,
          Real xmin, Real xmax,
          Real ymin, Real ymax,
          Real zmin, Real zmax,
          Real gamma,
          std::function<void(const Rvec&, const Rvec&, const Rvec&,
                             Vector3DR*, // density
                             Vector3DR*, // u
                             Vector3DR*, // v
                             Vector3DR*, // w
                             Vector3DR* /* pressure*/)> init_func)
  {
    Initialize(size_x, size_y, size_z,
               xmin, xmax, ymin, ymax, zmin, zmax,
               gamma, init_func);
  }

  void MainLoop(Int32 loop_max, Int32 out_step, const std::string &fname_base);
  void Update();
  void Output(const std::string &filename, Int32 step) const;
  void OutputSliceX(const std::string &filename, Int32 step) const;
  void OutputSliceY(const std::string &filename, Int32 step) const;
  void OutputSliceZ(const std::string &filename, Int32 step) const;
//  void OutputSiaX(const std::string &filename, Int32 step) const;
//  void OutputSiaY(const std::string &filename, Int32 step) const;
//  void OutputSiaZ(const std::string &filename, Int32 step) const;

  Fluid3D& Initialize(Int32 size_x, Int32 size_y, Int32 size_z,
                      Real xmin, Real xmax,
                      Real ymin, Real ymax,
                      Real zmin, Real zmax,
                      Real gamma,
                      std::function<void(const Rvec&, const Rvec&, const Rvec&,
                                         Vector3DR*,
                                         Vector3DR*,
                                         Vector3DR*,
                                         Vector3DR*,
                                         Vector3DR*)> init_func);

  void SetBCXLeft(std::function<void
                  (Int32 /*is*/, Int32 /*ie*/,
                   Int32 /*js*/, Int32 /*je*/,
                   Int32 /*ks*/, Int32 /*ke*/,
                   Vector3DR */*d_sia*/,
                   Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                   Vector3DR */*p_sia*/,
                   Vector3DR */*d_via*/,
                   Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                   Vector3DR */*p_via*/) noexcept> bc)
  { bc_x_left_ = bc; }

  void SetBCXRight(std::function<void
                   (Int32 /*is*/, Int32 /*ie*/,
                    Int32 /*js*/, Int32 /*je*/,
                    Int32 /*ks*/, Int32 /*ke*/,
                    Vector3DR */*d_sia*/,
                    Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                    Vector3DR */*p_sia*/,
                    Vector3DR */*d_via*/,
                    Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                    Vector3DR */*p_via*/) noexcept> bc)
  { bc_x_right_ = bc; }

  void SetBCYLeft(std::function<void
                  (Int32 /*is*/, Int32 /*ie*/,
                   Int32 /*js*/, Int32 /*je*/,
                   Int32 /*ks*/, Int32 /*ke*/,
                   Vector3DR */*d_sia*/,
                   Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                   Vector3DR */*p_sia*/,
                   Vector3DR */*d_via*/,
                   Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                   Vector3DR */*p_via*/) noexcept> bc)
  { bc_y_left_ = bc; }

  void SetBCYRight(std::function<void
                   (Int32 /*is*/, Int32 /*ie*/,
                    Int32 /*js*/, Int32 /*je*/,
                    Int32 /*ks*/, Int32 /*ke*/,
                    Vector3DR */*d_sia*/,
                    Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                    Vector3DR */*p_sia*/,
                    Vector3DR */*d_via*/,
                    Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                    Vector3DR */*p_via*/) noexcept> bc)
  { bc_y_right_ = bc; }

  void SetBCZLeft(std::function<void
                  (Int32 /*is*/, Int32 /*ie*/,
                   Int32 /*js*/, Int32 /*je*/,
                   Int32 /*ks*/, Int32 /*ke*/,
                   Vector3DR */*d_sia*/,
                   Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                   Vector3DR */*p_sia*/,
                   Vector3DR */*d_via*/,
                   Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                   Vector3DR */*p_via*/) noexcept> bc)
  { bc_z_left_ = bc; }

  void SetBCZRight(std::function<void
                   (Int32 /*is*/, Int32 /*ie*/,
                    Int32 /*js*/, Int32 /*je*/,
                    Int32 /*ks*/, Int32 /*ke*/,
                    Vector3DR */*d_sia*/,
                    Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                    Vector3DR */*p_sia*/,
                    Vector3DR */*d_via*/,
                    Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                    Vector3DR */*p_via*/) noexcept> bc)
  { bc_z_right_ = bc; }

private:
  Int32 nx_, ny_, nz_;
  Real xmin_, xmax_;
  Real ymin_, ymax_;
  Real zmin_, zmax_;
  Rvec x_;
  Rvec y_;
  Rvec z_;
  Real dt_;
  Real time_;

  Real gamma_{1.4};

  Vector3DR density_pv_;
  Vector3DR vx_pv_;
  Vector3DR vy_pv_;
  Vector3DR vz_pv_;
  Vector3DR pressure_pv_;

  Vector3DR density_lia_x_;
  Vector3DR vx_lia_x_;
  Vector3DR vy_lia_x_;
  Vector3DR vz_lia_x_;
  Vector3DR pressure_lia_x_;

  Vector3DR density_lia_y_;
  Vector3DR vx_lia_y_;
  Vector3DR vy_lia_y_;
  Vector3DR vz_lia_y_;
  Vector3DR pressure_lia_y_;

  Vector3DR density_lia_z_;
  Vector3DR vx_lia_z_;
  Vector3DR vy_lia_z_;
  Vector3DR vz_lia_z_;
  Vector3DR pressure_lia_z_;

  Vector3DR density_sia_x_;
  Vector3DR vx_sia_x_;
  Vector3DR vy_sia_x_;
  Vector3DR vz_sia_x_;
  Vector3DR pressure_sia_x_;

  Vector3DR density_sia_y_;
  Vector3DR vx_sia_y_;
  Vector3DR vy_sia_y_;
  Vector3DR vz_sia_y_;
  Vector3DR pressure_sia_y_;

  Vector3DR density_sia_z_;
  Vector3DR vx_sia_z_;
  Vector3DR vy_sia_z_;
  Vector3DR vz_sia_z_;
  Vector3DR pressure_sia_z_;

  Vector3DR density_via_;
  Vector3DR vx_via_;
  Vector3DR vy_via_;
  Vector3DR vz_via_;
  Vector3DR pressure_via_;

  Vector3DR density_pv_new_;
  Vector3DR vx_pv_new_;
  Vector3DR vy_pv_new_;
  Vector3DR vz_pv_new_;
  Vector3DR pressure_pv_new_;

  Vector3DR density_lia_x_new_;
  Vector3DR vx_lia_x_new_;
  Vector3DR vy_lia_x_new_;
  Vector3DR vz_lia_x_new_;
  Vector3DR pressure_lia_x_new_;

  Vector3DR density_lia_y_new_;
  Vector3DR vx_lia_y_new_;
  Vector3DR vy_lia_y_new_;
  Vector3DR vz_lia_y_new_;
  Vector3DR pressure_lia_y_new_;

  Vector3DR density_lia_z_new_;
  Vector3DR vx_lia_z_new_;
  Vector3DR vy_lia_z_new_;
  Vector3DR vz_lia_z_new_;
  Vector3DR pressure_lia_z_new_;

  Vector3DR density_sia_x_new_;
  Vector3DR vx_sia_x_new_;
  Vector3DR vy_sia_x_new_;
  Vector3DR vz_sia_x_new_;
  Vector3DR pressure_sia_x_new_;

  Vector3DR density_sia_y_new_;
  Vector3DR vx_sia_y_new_;
  Vector3DR vy_sia_y_new_;
  Vector3DR vz_sia_y_new_;
  Vector3DR pressure_sia_y_new_;

  Vector3DR density_sia_z_new_;
  Vector3DR vx_sia_z_new_;
  Vector3DR vy_sia_z_new_;
  Vector3DR vz_sia_z_new_;
  Vector3DR pressure_sia_z_new_;

  Vector3DR density_via_new_;
  Vector3DR vx_via_new_;
  Vector3DR vy_via_new_;
  Vector3DR vz_via_new_;
  Vector3DR pressure_via_new_;

  Vector3DR dflux_;
  Vector3DR lflux_;
  Vector3DR mflux_;
  Vector3DR nflux_;
  Vector3DR eflux_;

  std::function<void(Int32 /*is*/, Int32 /*ie*/,
                     Int32 /*js*/, Int32 /*je*/,
                     Int32 /*ks*/, Int32 /*ke*/,
                     Vector3DR */*d_sia*/,
                     Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                     Vector3DR */*p_sia*/,
                     Vector3DR */*d_via*/,
                     Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                     Vector3DR */*p_via*/) noexcept> bc_x_left_;

  std::function<void(Int32 /*is*/, Int32 /*ie*/,
                     Int32 /*js*/, Int32 /*je*/,
                     Int32 /*ks*/, Int32 /*ke*/,
                     Vector3DR */*d_sia*/,
                     Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                     Vector3DR */*p_sia*/,
                     Vector3DR */*d_via*/,
                     Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                     Vector3DR */*p_via*/) noexcept> bc_x_right_;

  std::function<void(Int32 /*is*/, Int32 /*ie*/,
                     Int32 /*js*/, Int32 /*je*/,
                     Int32 /*ks*/, Int32 /*ke*/,
                     Vector3DR */*d_sia*/,
                     Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                     Vector3DR */*p_sia*/,
                     Vector3DR */*d_via*/,
                     Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                     Vector3DR */*p_via*/) noexcept> bc_y_left_;

  std::function<void(Int32 /*is*/, Int32 /*ie*/,
                     Int32 /*js*/, Int32 /*je*/,
                     Int32 /*ks*/, Int32 /*ke*/,
                     Vector3DR */*d_sia*/,
                     Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                     Vector3DR */*p_sia*/,
                     Vector3DR */*d_via*/,
                     Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                     Vector3DR */*p_via*/) noexcept> bc_y_right_;

  std::function<void(Int32 /*is*/, Int32 /*ie*/,
                     Int32 /*js*/, Int32 /*je*/,
                     Int32 /*ks*/, Int32 /*ke*/,
                     Vector3DR */*d_sia*/,
                     Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                     Vector3DR */*p_sia*/,
                     Vector3DR */*d_via*/,
                     Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                     Vector3DR */*p_via*/) noexcept> bc_z_left_;

  std::function<void(Int32 /*is*/, Int32 /*ie*/,
                     Int32 /*js*/, Int32 /*je*/,
                     Int32 /*ks*/, Int32 /*ke*/,
                     Vector3DR */*d_sia*/,
                     Vector3DR */*u_sia*/, Vector3DR */*v_sia*/, Vector3DR */*w_sia*/,
                     Vector3DR */*p_sia*/,
                     Vector3DR */*d_via*/,
                     Vector3DR */*u_via*/, Vector3DR */*v_via*/, Vector3DR */*w_via*/,
                     Vector3DR */*p_via*/) noexcept> bc_z_right_;

  void NewDt() noexcept;

  void BoundaryConditions() noexcept;
  void BoundaryConditionsX() noexcept;
  void BoundaryConditionsY() noexcept;
  void BoundaryConditionsZ() noexcept;

  void Advection();
  void AdvectionX(Int32 ny, Int32 nz,
                  const Vector3DR& d_p, const Vector3DR& u_p,
                  const Vector3DR& v_p, const Vector3DR& w_p,
                  const Vector3DR& p_p, const Vector3DR& d_v,
                  const Vector3DR& u_v, const Vector3DR& v_v,
                  const Vector3DR& w_v, const Vector3DR& p_v,
                  Vector3DR* d_pn, Vector3DR* u_pn,
                  Vector3DR* v_pn, Vector3DR* w_pn,
                  Vector3DR* p_pn, Vector3DR* d_vn,
                  Vector3DR* u_vn, Vector3DR* v_vn,
                  Vector3DR* w_vn,Vector3DR* p_vn) noexcept;

  void AdvectionY(Int32 nz, Int32 nx,
                  const Vector3DR& d_p, const Vector3DR& u_p,
                  const Vector3DR& v_p, const Vector3DR& w_p,
                  const Vector3DR& p_p, const Vector3DR& d_v,
                  const Vector3DR& u_v, const Vector3DR& v_v,
                  const Vector3DR& w_v, const Vector3DR& p_v,
                  Vector3DR* d_pn, Vector3DR* u_pn,
                  Vector3DR* v_pn, Vector3DR* w_pn,
                  Vector3DR* p_pn, Vector3DR* d_vn,
                  Vector3DR* u_vn, Vector3DR* v_vn,
                  Vector3DR* w_vn,Vector3DR* p_vn) noexcept;

  void AdvectionZ(Int32 nx, Int32 ny,
                  const Vector3DR& d_p, const Vector3DR& u_p,
                  const Vector3DR& v_p, const Vector3DR& w_p,
                  const Vector3DR& p_p, const Vector3DR& d_v,
                  const Vector3DR& u_v, const Vector3DR& v_v,
                  const Vector3DR& w_v, const Vector3DR& p_v,
                  Vector3DR* d_pn, Vector3DR* u_pn,
                  Vector3DR* v_pn, Vector3DR* w_pn,
                  Vector3DR* p_pn, Vector3DR* d_vn,
                  Vector3DR* u_vn, Vector3DR* v_vn,
                  Vector3DR* w_vn,Vector3DR* p_vn) noexcept;

  void AdvCipRK3X(Int32 l, Int32 i, Int32 j, Int32 k,
                  Real dens[4], Real velx[4], Real vely[4], Real velz[4], Real pres[4],
                  const Vector3DR& d_p, const Vector3DR& u_p,
                  const Vector3DR& v_p, const Vector3DR& w_p,
                  const Vector3DR& p_p, const Vector3DR& d_v,
                  const Vector3DR& u_v, const Vector3DR& v_v,
                  const Vector3DR& w_v, const Vector3DR& p_v) noexcept;

  void AdvCipRK3Y(Int32 l, Int32 i, Int32 j, Int32 k,
                  Real dens[4], Real velx[4], Real vely[4], Real velz[4], Real pres[4],
                  const Vector3DR& d_p, const Vector3DR& u_p,
                  const Vector3DR& v_p, const Vector3DR& w_p,
                  const Vector3DR& p_p, const Vector3DR& d_v,
                  const Vector3DR& u_v, const Vector3DR& v_v,
                  const Vector3DR& w_v, const Vector3DR& p_v) noexcept;

  void AdvCipRK3Z(Int32 l, Int32 i, Int32 j, Int32 k,
                  Real dens[4], Real velx[4], Real vely[4], Real velz[4], Real pres[4],
                  const Vector3DR& d_p, const Vector3DR& u_p,
                  const Vector3DR& v_p, const Vector3DR& w_p,
                  const Vector3DR& p_p, const Vector3DR& d_v,
                  const Vector3DR& u_v, const Vector3DR& v_v,
                  const Vector3DR& w_v, const Vector3DR& p_v) noexcept;

  Real CipCsl3X(Int32 i, Int32 j, Int32 k, Real v,
                const Vector3DR& sia, const Vector3DR& via) noexcept;
  Real CipCsl3Y(Int32 i, Int32 j, Int32 k, Real v,
                const Vector3DR& sia, const Vector3DR& via) noexcept;
  Real CipCsl3Z(Int32 i, Int32 j, Int32 k, Real v,
                const Vector3DR& sia, const Vector3DR& via) noexcept;

  void UpdateValue();
};

#endif // FLUID3D_H
