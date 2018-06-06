#include <string>

#include "global_definitions.h"
#include "Fluid3D.h"

void Fluid3D::MainLoop(Int32 loop_max, Int32 out_step,
                       const std::string &fname_base)
{
  Int32 output = 0;

  for(auto step = 0; step < loop_max; ++step){
    if(step % out_step == 0){
      Output(fname_base, output);
      OutputSliceX(fname_base, output);
      OutputSliceY(fname_base, output);
      OutputSliceZ(fname_base, output);
      ++output;
    }
    Update();
  }
  Output(fname_base, output);
  OutputSliceX(fname_base, output);
  OutputSliceY(fname_base, output);
  OutputSliceZ(fname_base, output);
}

void Fluid3D::Update()
{
  Advection();
  BoundaryConditions();
  time_ += dt_;
  NewDt();
}

void Fluid3D::BoundaryConditions() noexcept
{
  BoundaryConditionsX();
  BoundaryConditionsY();
  BoundaryConditionsZ();
}

void Fluid3D::UpdateValue()
{
  // PV
  swap(density_pv_, density_pv_new_);
  swap(vx_pv_, vx_pv_new_);
  swap(vy_pv_, vy_pv_new_);
  swap(vz_pv_, vz_pv_new_);
  swap(pressure_pv_, pressure_pv_new_);

  // LIA x
  swap(density_lia_x_, density_lia_x_new_);
  swap(vx_lia_x_, vx_lia_x_new_);
  swap(vy_lia_x_, vy_lia_x_new_);
  swap(vz_lia_x_, vz_lia_x_new_);
  swap(pressure_lia_x_, pressure_lia_x_new_);

  // LIA y
  swap(density_lia_y_, density_lia_y_new_);
  swap(vx_lia_y_, vx_lia_y_new_);
  swap(vy_lia_y_, vy_lia_y_new_);
  swap(vz_lia_y_, vz_lia_y_new_);
  swap(pressure_lia_y_, pressure_lia_y_new_);

  // LIA z
  swap(density_lia_z_, density_lia_z_new_);
  swap(vx_lia_z_, vx_lia_z_new_);
  swap(vy_lia_z_, vy_lia_z_new_);
  swap(vz_lia_z_, vz_lia_z_new_);
  swap(pressure_lia_z_, pressure_lia_z_new_);

  // SIA x
  swap(density_sia_x_, density_sia_x_new_);
  swap(vx_sia_x_, vx_sia_x_new_);
  swap(vy_sia_x_, vy_sia_x_new_);
  swap(vz_sia_x_, vz_sia_x_new_);
  swap(pressure_sia_x_, pressure_sia_x_new_);

  // SIA y
  swap(density_sia_y_, density_sia_y_new_);
  swap(vx_sia_y_, vx_sia_y_new_);
  swap(vy_sia_y_, vy_sia_y_new_);
  swap(vz_sia_y_, vz_sia_y_new_);
  swap(pressure_sia_y_, pressure_sia_y_new_);

  // SIA z
  swap(density_sia_z_, density_sia_z_new_);
  swap(vx_sia_z_, vx_sia_z_new_);
  swap(vy_sia_z_, vy_sia_z_new_);
  swap(vz_sia_z_, vz_sia_z_new_);
  swap(pressure_sia_z_, pressure_sia_z_new_);

  // VIA
  swap(density_via_, density_via_new_);
  swap(vx_via_, vx_via_new_);
  swap(vy_via_, vy_via_new_);
  swap(vz_via_, vz_via_new_);
  swap(pressure_via_, pressure_via_new_);
}
