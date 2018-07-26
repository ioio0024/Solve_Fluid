#ifndef FLUID3D_H
#define FLUID3D_H

#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>

#include "global_definitions.h"
#include "FluidQuantites.h"
#include "Vector3D.h"

template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
class Fluid3D
{
public:
  Fluid3D() = default;
  Fluid3D(Int32 size_x, Int32 size_y, Int32 size_z,
          Real xmin, Real xmax,
          Real ymin, Real ymax,
          Real zmin, Real zmax,
          Real gamma,
          std::function<void(const Vector3DR&, const Vector3DR&, const Vector3DR&,
                             FluidQuantity3DArray*)> init_func)
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

  void Initialize(Int32 size_x, Int32 size_y, Int32 size_z,
                  Real xmin, Real xmax,
                  Real ymin, Real ymax,
                  Real zmin, Real zmax,
                  Real gamma,
                  std::function<void(const Vector3DR&, const Vector3DR&, const Vector3DR&,
                                     FluidQuantity3DArray*)> init_func);

private:
  Int32 nx_, ny_, nz_;
  Real xmin_, xmax_;
  Real ymin_, ymax_;
  Real zmin_, zmax_;
  Vector3DR x_;
  Vector3DR y_;
  Vector3DR z_;
  Real dt_;
  Real time_;

  Real gamma_{1.4};

  FluidQuantity3DArray pv_;

  FluidQuantity3DArray lia_x_;
  FluidQuantity3DArray lia_y_;
  FluidQuantity3DArray lia_z_;

  FluidQuantity3DArray sia_x_;
  FluidQuantity3DArray sia_y_;
  FluidQuantity3DArray sia_z_;

  FluidQuantity3DArray via_;

  FluidQuantity3DArray pv_new_;

  FluidQuantity3DArray lia_x_new_;
  FluidQuantity3DArray lia_y_new_;
  FluidQuantity3DArray lia_z_new_;

  FluidQuantity3DArray sia_x_new_;
  FluidQuantity3DArray sia_y_new_;
  FluidQuantity3DArray sia_z_new_;

  FluidQuantity3DArray via_new_;

  FluidQuantity3DArray flux_;

  BCxl bc_x_left_;
  BCxr bc_x_right_;

  BCyl bc_y_left_;
  BCyr bc_y_right_;

  BCzl bc_z_left_;
  BCzr bc_z_right_;

  CIPx cipx_;
  CIPy cipy_;
  CIPz cipz_;

  void NewDt() noexcept;

  void BoundaryConditions() noexcept;
  void BoundaryConditionsX() noexcept;
  void BoundaryConditionsY() noexcept;
  void BoundaryConditionsZ() noexcept;

  void Advection();
  void AdvectionX(Int32 ny, Int32 nz,
                  const FluidQuantity3DArray &sia,
                  const FluidQuantity3DArray &via,
                  FluidQuantity3DArray *sia_new,
                  FluidQuantity3DArray *via_new
                  ) noexcept;

  void AdvectionY(Int32 nz, Int32 nx,
                  const FluidQuantity3DArray &sia,
                  const FluidQuantity3DArray &via,
                  FluidQuantity3DArray *sia_new,
                  FluidQuantity3DArray *via_new
                  ) noexcept;

  void AdvectionZ(Int32 nx, Int32 ny,
                  const FluidQuantity3DArray &sia,
                  const FluidQuantity3DArray &via,
                  FluidQuantity3DArray *sia_new,
                  FluidQuantity3DArray *via_new
                  ) noexcept;

  void AdvCipRK3X(Int32 l, Int32 i, Int32 j, Int32 k,
                  FluidQuantity3D rk[4],
                  const FluidQuantity3DArray &sia,
                  const FluidQuantity3DArray &via
                  ) noexcept;

  void AdvCipRK3Y(Int32 l, Int32 i, Int32 j, Int32 k,
                  FluidQuantity3D rk[4],
                  const FluidQuantity3DArray &sia,
                  const FluidQuantity3DArray &via
                  ) noexcept;

  void AdvCipRK3Z(Int32 l, Int32 i, Int32 j, Int32 k,
                  FluidQuantity3D rk[4],
                  const FluidQuantity3DArray &sia,
                  const FluidQuantity3DArray &via
                  ) noexcept;

  void UpdateValue();
};

#include "Advection.h"
#include "BoundaryConditions.h"
#include "Fluid3DMain.h"
#include "Initialize.h"
#include "NewDt.h"
#include "Output.h"

#endif // FLUID3D_H
