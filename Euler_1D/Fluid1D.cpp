#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "global_definitions.h"
#include "Fluid1D.h"

Fluid1D::Fluid1D(Int32 size, Real xmin, Real xmax, Real gamma,
                 std::function<void (const Rvec &,
                                     Rvec *,
                                     Rvec *,
                                     Rvec *)> init_func)
  : size_{size + 2 * GHOST}, xmin_{xmin}, xmax_{xmax},
    x_(size_ + 1), dt_{1.0e3}, time_{0.0}, gamma_{gamma},
    density_prim_sia_(size_ + 1), density_prim_sia_new_(size_ + 1),
    velocity_prim_sia_(size_ + 1), velocity_prim_sia_new_(size_ + 1),
    pressure_prim_sia_(size_ + 1), pressure_prim_sia_new_(size_ + 1),
    density_prim_via_(size_), density_prim_via_new_(size_),
    velocity_prim_via_(size_), velocity_prim_via_new_(size_),
    pressure_prim_via_(size_), pressure_prim_via_new_(size_),
    density_cons_via_(size_), density_flux_(size_ + 1),
    momentum_cons_via_(size_), momentum_flux_(size_ + 1),
    energy_cons_via_(size_), energy_flux_(size_ + 1)
{
  auto dx = (xmax_ - xmin_) / size;

  for(auto i = 0; i < size_ + 1; ++i){
    x_[i] = xmin_ + dx * (static_cast<double>(i) - static_cast<double>(GHOST));
  }

  init_func(x_,
            &density_prim_sia_,
            &velocity_prim_sia_,
            &pressure_prim_sia_);

  for(auto i = 0; i < size_; ++i){
    density_prim_via_[i] = 0.5 * (density_prim_sia_[i + 1]
                                + density_prim_sia_[i]);
    velocity_prim_via_[i] = 0.5 * (velocity_prim_sia_[i + 1]
                                 +velocity_prim_sia_[i]);
    pressure_prim_via_[i] = 0.5 * (pressure_prim_sia_[i + 1]
                                 + pressure_prim_sia_[i]);
  }

  for(auto i = 0; i < size_; ++i) {
    density_cons_via_[i] = density_prim_via_[i];
    momentum_cons_via_[i] = density_prim_via_[i] * velocity_prim_via_[i];
    energy_cons_via_[i] = pressure_prim_via_[i] / (gamma_ - 1.0)
                        + 0.5 * density_prim_via_[i]
                        * velocity_prim_via_[i] * velocity_prim_via_[i];
  }

  BoundaryConditions();

  NewDt();
  dt_ *= 0.01;
}

Fluid1D& Fluid1D::Initialize(Int32 size, Real xmin, Real xmax, Real gamma,
                         std::function<void (const Rvec &,
                                             Rvec *,
                                             Rvec *,
                                             Rvec *)> init_func)
{
  Fluid1D tmp(size, xmin, xmax, gamma, init_func);
  *this = tmp;
  return *this;
}

void Fluid1D::MainLoop(Int32 loop_max, Int32 out_step,
                       const std::string &fname_base)
{
    Int32 output = 0;

    for(auto step = 0; step < loop_max; ++step){
        if(step % out_step == 0){
            Output(fname_base, output++);
        }
        Update();
    }
    Output(fname_base, output);
}

void Fluid1D::Update() noexcept
{
  Advection();
  UpdateValue();
  BoundaryConditions();
  time_ += dt_;
  NewDt();
}

void Fluid1D::NewDt() noexcept
{
  Real vmax = 1.0e-10;
  for(auto i = GHOST; i < size_ - GHOST; ++i){
    auto v = velocity_prim_via_[i];
    auto cs = std::sqrt(gamma_ * pressure_prim_via_[i] / density_prim_via_[i]);
    auto dx = x_[i + 1] - x_[i];

    vmax = std::fmax(vmax, std::fabs(v));
    vmax = std::fmax(vmax, std::fabs(v + cs));
    vmax = std::fmax(vmax, std::fabs(v - cs));
    dt_ = std::fmin(dt_, CFL * dx / (vmax + EPS));
  }
}

void Fluid1D::BoundaryConditions() noexcept
{
  auto &dps = density_prim_sia_;
  auto &vps = velocity_prim_sia_;
  auto &pps = pressure_prim_sia_;
  auto &dpv = density_prim_via_;
  auto &vpv = velocity_prim_via_;
  auto &ppv = pressure_prim_via_;
  for(auto i = 0; i <= GHOST; ++i){
    dps[i] = dpv[GHOST];
    vps[i] = 0.0;
    pps[i] = ppv[GHOST];

    dps[size_ - GHOST + i] = dpv[size_ - 1 - GHOST];
    vps[size_ - GHOST + i] = 0.0;
    pps[size_ - GHOST + i] = ppv[size_ - 1 - GHOST];
  }

  auto &dcv = density_cons_via_;
  auto &mcv = momentum_cons_via_;
  auto &ecv = energy_cons_via_;
  for(auto i = 0; i < GHOST; ++i){
    dpv[i] = dpv[GHOST];
    vpv[i] = -vpv[GHOST];
    ppv[i] = ppv[GHOST];

    dcv[i] = dcv[GHOST];
    mcv[i] = -mcv[GHOST];
    ecv[i] = ecv[GHOST];

    dpv[size_ - GHOST + i] = dpv[size_ - GHOST - 1];
    vpv[size_ - GHOST + i] = -vpv[size_ - GHOST - 1];
    ppv[size_ - GHOST + i] = ppv[size_ - GHOST - 1];

    dcv[size_ - GHOST + i] = dcv[size_ - GHOST - 1];
    mcv[size_ - GHOST + i] = -mcv[size_ - GHOST - 1];
    ecv[size_ - GHOST + i] = ecv[size_ - GHOST - 1];
  }
}

void Fluid1D::UpdateValue() noexcept
{
  auto &dps = density_prim_sia_;
  auto &dpsn = density_prim_sia_new_;
  auto &vps = velocity_prim_sia_;
  auto &vpsn = velocity_prim_sia_new_;
  auto &pps = pressure_prim_sia_;
  auto &ppsn = pressure_prim_sia_new_;
  auto &dpv = density_prim_via_;
  auto &dpvn = density_prim_via_new_;
  auto &vpv = velocity_prim_via_;
  auto &vpvn = velocity_prim_via_new_;
  auto &ppv = pressure_prim_via_;
  auto &ppvn = pressure_prim_via_new_;

  for(auto i = 0; i < size_ + 1; ++i){
    dps[i] = dpsn[i];
    vps[i] = vpsn[i];
    pps[i] = ppsn[i];
  }
  for(auto i = 0; i < size_; ++i){
    dpv[i] = dpvn[i];
    vpv[i] = vpvn[i];
    ppv[i] = ppvn[i];
  }
}

void Fluid1D::Output(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << '-' << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto i = GHOST; i < size_ - GHOST; ++i){
    file << (x_[i] + x_[i + 1]) * 0.5 << ' ' << density_cons_via_[i] << ' '
         << momentum_cons_via_[i]  << ' ' << energy_cons_via_[i] << ' '
         << density_prim_via_[i] << ' ' << velocity_prim_via_[i] << ' '
         << pressure_prim_via_[i] << '\n';
  }

  file.close ();
}
