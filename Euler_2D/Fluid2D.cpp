#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "global_definitions.h"
#include "Fluid2D.h"

Fluid2D& Fluid2D::Initialize(Int32 size_x, Int32 size_y,
                    Real xmin, Real xmax, Real ymin, Real ymax, Real gamma,
                    std::function<void (const Rvec &, // x
                                        const Rvec &, // y
                                        Vector2D *,   // density
                                        Vector2D *,   // vx
                                        Vector2D *,   // vy
                                        Vector2D *    // pressure
                                        )> init_func)
{
  size_x_ = size_x + 2 * GHOST;
  size_y_ = size_y + 2 * GHOST;
  xmin_ = xmin;
  xmax_ = xmax;
  ymin_ = ymin;
  ymax_ = ymax;

  x_.resize(size_x_ + 1);
  y_.resize(size_y_ + 1);

  dt_ = 1.0e3;
  time_ = 0.0;
  gamma_ = gamma;

  density_prim_pv_.resize(size_x_ + 1, size_y_ + 1);
  vx_prim_pv_.resize(size_x_ + 1, size_y_ + 1);
  vy_prim_pv_.resize(size_x_ + 1, size_y_ + 1);
  pressure_prim_pv_.resize(size_x_ + 1, size_y_ + 1);

  density_prim_sia_x_.resize(size_x_ + 1, size_y_);
  vx_prim_sia_x_.resize(size_x_ + 1, size_y_);
  vy_prim_sia_x_.resize(size_x_ + 1, size_y_);
  pressure_prim_sia_x_.resize(size_x_ + 1, size_y_);

  density_prim_sia_y_.resize(size_x_, size_y_ + 1);
  vx_prim_sia_y_.resize(size_x_, size_y_ + 1);
  vy_prim_sia_y_.resize(size_x_, size_y_ + 1);
  pressure_prim_sia_y_.resize(size_x_, size_y_ + 1);

  density_prim_via_.resize(size_x_, size_y_);
  vx_prim_via_.resize(size_x_, size_y_);
  vy_prim_via_.resize(size_x_, size_y_);
  pressure_prim_via_.resize(size_x_, size_y_);

  density_prim_pv_new_.resize(size_x_ + 1, size_y_ + 1);
  vx_prim_pv_new_.resize(size_x_ + 1, size_y_ + 1);
  vy_prim_pv_new_.resize(size_x_ + 1, size_y_ + 1);
  pressure_prim_pv_new_.resize(size_x_ + 1, size_y_ + 1);

  density_prim_sia_x_new_.resize(size_x_ + 1, size_y_);
  vx_prim_sia_x_new_.resize(size_x_ + 1, size_y_);
  vy_prim_sia_x_new_.resize(size_x_ + 1, size_y_);
  pressure_prim_sia_x_new_.resize(size_x_ + 1, size_y_);

  density_prim_sia_y_new_.resize(size_x_, size_y_ + 1);
  vx_prim_sia_y_new_.resize(size_x_, size_y_ + 1);
  vy_prim_sia_y_new_.resize(size_x_, size_y_ + 1);
  pressure_prim_sia_y_new_.resize(size_x_, size_y_ + 1);

  density_prim_via_new_.resize(size_x_, size_y_);
  vx_prim_via_new_.resize(size_x_, size_y_);
  vy_prim_via_new_.resize(size_x_, size_y_);
  pressure_prim_via_new_.resize(size_x_, size_y_);

  dflux_x_.resize(size_x_ + 1, size_y_ + 1);
  lflux_x_.resize(size_x_ + 1, size_y_ + 1);
  mflux_x_.resize(size_x_ + 1, size_y_ + 1);
  eflux_x_.resize(size_x_ + 1, size_y_ + 1);

  dflux_y_.resize(size_x_ + 1, size_y_ + 1);
  lflux_y_.resize(size_x_ + 1, size_y_ + 1);
  mflux_y_.resize(size_x_ + 1, size_y_ + 1);
  eflux_y_.resize(size_x_ + 1, size_y_ + 1);

  auto dx = (xmax_ - xmin_) / size_x;
  auto dy = (ymax_ - ymin_) / size_y;

  for(auto i = 0; i < size_x_ + 1; ++i){
    x_[i] = xmin_ + dx * (static_cast<Real>(i) - static_cast<Real>(GHOST));
  }
  for(auto j = 0; j < size_y_ + 1; ++j){
    y_[j] = ymin_ + dy * (static_cast<Real>(j) - static_cast<Real>(GHOST));
  }

  init_func(x_, y_,
            &density_prim_pv_,
            &vx_prim_pv_,
            &vy_prim_pv_,
            &pressure_prim_pv_);

  // sia_x
  for(auto j = 0; j < size_y_; ++j){
    for(auto i = 0; i < size_x_ + 1; ++i){
      density_prim_sia_x_(i, j) = 0.5 * (
            density_prim_pv_(i, j) + density_prim_pv_(i, j + 1));
      vx_prim_sia_x_(i, j) = 0.5 * (
            vx_prim_pv_(i, j) + vx_prim_pv_(i, j + 1));
      vy_prim_sia_x_(i, j) = 0.5 * (
            vy_prim_pv_(i, j) + vy_prim_pv_(i, j + 1));
      pressure_prim_sia_x_(i, j) = 0.5 * (
            pressure_prim_pv_(i, j) + pressure_prim_pv_(i, j + 1));
    }
  }

  // sia_y
  for(auto j = 0; j < size_y_ + 1; ++j){
    for(auto i = 0; i < size_x_; ++i){
      density_prim_sia_y_(i, j) = 0.5 * (
            density_prim_pv_(i, j) + density_prim_pv_(i + 1, j));
      vx_prim_sia_y_(i, j) = 0.5 * (
            vx_prim_pv_(i, j) + vx_prim_pv_(i + 1, j));
      vy_prim_sia_y_(i, j) = 0.5 * (
            vy_prim_pv_(i, j) + vy_prim_pv_(i + 1, j));
      pressure_prim_sia_y_(i, j) = 0.5 * (
            pressure_prim_pv_(i, j) + pressure_prim_pv_(i + 1, j));
    }
  }

  // via
  for(auto j = 0; j < size_y_; ++j){
    for(auto i = 0; i < size_x_; ++i){
      density_prim_via_(i, j) = 0.25 * (
            density_prim_pv_(i, j) + density_prim_pv_(i + 1, j) +
            density_prim_pv_(i, j + 1) + density_prim_pv_(i + 1, j + 1));
      vx_prim_via_(i, j) = 0.25 * (
            vx_prim_pv_(i, j) + vx_prim_pv_(i + 1, j) +
            vx_prim_pv_(i, j + 1) + vx_prim_pv_(i + 1, j + 1));
      vy_prim_via_(i, j) = 0.25 * (
            vy_prim_pv_(i, j) + vy_prim_pv_(i + 1, j) +
            vy_prim_pv_(i, j + 1) + vy_prim_pv_(i + 1, j + 1));
      pressure_prim_via_(i, j) = 0.25 * (
            pressure_prim_pv_(i, j) + pressure_prim_pv_(i + 1, j) +
            pressure_prim_pv_(i, j + 1) + pressure_prim_pv_(i + 1, j + 1));
    }
  }

  BoundaryConditions();

  NewDt();
  dt_ *= 0.01;

  return *this;
}

void Fluid2D::MainLoop(Int32 loop_max, Int32 out_step,
                       const std::string &fname_base)
{
    Int32 output = 0;

    for(auto step = 0; step < loop_max; ++step){
        if(step % out_step == 0){
            Output(fname_base, output);
#ifndef NDEBUG
            OutputSiaX(fname_base, output);
            OutputSiaY(fname_base, output);
#endif
            ++output;
        }
        Update();
    }
    Output(fname_base, output);
#ifndef NDEBUG
    OutputSiaX(fname_base, output);
    OutputSiaY(fname_base, output);
#endif
}

void Fluid2D::Update() noexcept
{
  Advection();
  BoundaryConditions();
  time_ += dt_;
  NewDt();
}

void Fluid2D::NewDt() noexcept
{
  for(auto j = GHOST; j < size_y_ - GHOST + 1; ++j){
    Real vy, cs;
    auto dy = y_[j + 1] - y_[j];
    for(auto i = GHOST; i < size_x_ - GHOST + 1; ++i){
      auto vx = vx_prim_pv_(i, j);
      vy = vy_prim_pv_(i, j);
      cs = std::sqrt(gamma_
              * pressure_prim_pv_(i, j) / density_prim_pv_(i, j));
      auto dx = x_[i + 1] - x_[i];

      Real vmax = 1.0e-10;
      vmax = std::fmax(vmax, std::fabs(vx));
      vmax = std::fmax(vmax, std::fabs(vx + cs));
      vmax = std::fmax(vmax, std::fabs(vx - cs));
      dt_ = std::fmin(dt_, CFL * dx / (vmax + EPS));

      vmax = 1.0e-10;
      vmax = std::fmax(vmax, std::fabs(vy));
      vmax = std::fmax(vmax, std::fabs(vy + cs));
      vmax = std::fmax(vmax, std::fabs(vy - cs));
      dt_ = std::fmin(dt_, CFL * dy / (vmax + EPS));
    }
  }
/*
  for(auto j = GHOST; j < size_y_ - GHOST; ++j){
    Real vy, cs;
    auto dy = y_[j + 1] - y_[j];
    for(auto i = GHOST; i < size_x_ - GHOST; ++i){
      auto vx = vx_prim_via_(i, j);
      vy = vy_prim_via_(i, j);
      cs = std::sqrt(gamma_
              * pressure_prim_via_(i, j) / density_prim_via_(i, j));
      auto dx = x_[i + 1] - x_[i];

      Real vmax = 1.0e-10;
      vmax = std::fmax(vmax, std::fabs(vx));
      vmax = std::fmax(vmax, std::fabs(vx + cs));
      vmax = std::fmax(vmax, std::fabs(vx - cs));
      dt_ = std::fmin(dt_, CFL * dx / (vmax + EPS));

      vmax = std::fmax(vmax, std::fabs(vy));
      vmax = std::fmax(vmax, std::fabs(vy + cs));
      vmax = std::fmax(vmax, std::fabs(vy - cs));
      dt_ = std::fmin(dt_, CFL * dy / (vmax + EPS));
    }
  }
*/
}

void Fluid2D::BoundaryConditions() noexcept
{
  BoundaryConditionsX();
  BoundaryConditionsY();
}

void Fluid2D::BoundaryConditionsX() noexcept
{
  for(auto j = 0; j < size_y_ + 1; ++j){
    for(auto i = 0; i < GHOST + 1; ++i){
      // PV
      density_prim_pv_(i, j) = density_prim_pv_(GHOST, j);
      vx_prim_pv_(i, j) = 0.0;
      vy_prim_pv_(i, j) = vy_prim_pv_(GHOST, j);
      pressure_prim_pv_(i, j) = pressure_prim_pv_(GHOST, j);

      density_prim_pv_(size_x_ - GHOST + i, j)
                      = density_prim_pv_(size_x_ - GHOST, j);
      vx_prim_pv_(size_x_ - GHOST + i, j) = 0.0;
      vy_prim_pv_(size_x_ - GHOST + i, j)
                 = vy_prim_pv_(size_x_ - GHOST, j);
      pressure_prim_pv_(size_x_ - GHOST + i, j)
                       = pressure_prim_pv_(size_x_ - GHOST, j);
    }
  }

  // sia_x
  for(auto j = 0; j < size_y_; ++j){
    for(auto i = 0; i < GHOST + 1; ++i){
      density_prim_sia_x_(i, j) = density_prim_sia_x_(GHOST, j);
      vx_prim_sia_x_(i, j) = 0.0;
      vy_prim_sia_x_(i, j) = vy_prim_sia_x_(GHOST, j);
      pressure_prim_sia_x_(i, j) = pressure_prim_sia_x_(GHOST, j);

      density_prim_sia_x_(size_x_ - GHOST + i, j) =
          density_prim_sia_x_(size_x_ - GHOST, j);
      vx_prim_sia_x_(size_x_ - GHOST + i, j) = 0.0;
      vy_prim_sia_x_(size_x_ - GHOST + i, j) =
          vy_prim_sia_x_(size_x_ - GHOST, j);
      pressure_prim_sia_x_(size_x_ - GHOST + i, j) =
          pressure_prim_sia_x_(size_x_ - GHOST, j);
    }
  }

  // sia_y
  for(auto j = 0; j < size_y_ + 1; ++j){
    for(auto i = 0; i < GHOST; ++i){
      density_prim_sia_y_(i, j) = density_prim_sia_y_(GHOST, j);
      vx_prim_sia_y_(i, j) = -vx_prim_sia_y_(GHOST, j);
      vy_prim_sia_y_(i, j) = vy_prim_sia_y_(GHOST, j);
      pressure_prim_sia_y_(i, j) = pressure_prim_sia_y_(GHOST, j);

      density_prim_sia_y_(size_x_ - GHOST + i, j) =
          density_prim_sia_y_(size_x_ - GHOST - 1, j);
      vx_prim_sia_y_(size_x_ - GHOST + i, j) =
          -vx_prim_sia_y_(size_x_ - GHOST - 1, j);
      vy_prim_sia_y_(size_x_ - GHOST + i, j) =
          vy_prim_sia_y_(size_x_ - GHOST - 1, j);
      pressure_prim_sia_y_(size_x_ - GHOST + i, j) =
          pressure_prim_sia_y_(size_x_ - GHOST - 1, j);
    }
  }

  // via
  for(auto j = 0; j < size_y_; ++j){
    for(auto i = 0; i < GHOST; ++i){
      density_prim_via_(i, j) = density_prim_via_(GHOST, j);
      vx_prim_via_(i, j) = -vx_prim_via_(GHOST, j);
      vy_prim_via_(i, j) = vy_prim_via_(GHOST, j);
      pressure_prim_via_(i, j) = pressure_prim_via_(GHOST, j);

      density_prim_via_(size_x_ - GHOST + i, j) =
          density_prim_via_(size_x_ - GHOST - 1, j);
      vx_prim_via_(size_x_ - GHOST + i, j) =
          -vx_prim_via_(size_x_ - GHOST - 1, j);
      vy_prim_via_(size_x_ - GHOST + i, j) =
          vy_prim_via_(size_x_ - GHOST - 1, j);
      pressure_prim_via_(size_x_ - GHOST + i, j) =
          pressure_prim_via_(size_x_ - GHOST - 1, j);
    }
  }
}

void Fluid2D::BoundaryConditionsY() noexcept
{
  for(auto j = 0; j < GHOST + 1; ++j){
    for(auto i = 0; i < size_x_ + 1; ++i){
      // PV
      density_prim_pv_(i, j) = density_prim_pv_(i, GHOST);
      vx_prim_pv_(i, j) = vx_prim_pv_(i, GHOST);
      vy_prim_pv_(i, j) = 0.0;
      pressure_prim_pv_(i, j) = pressure_prim_pv_(i, GHOST);

      density_prim_pv_(i, size_y_ - j)
                      = density_prim_pv_(i, size_y_ - GHOST);
      vx_prim_pv_(i, size_y_ - j) =
          vx_prim_pv_(i, size_y_ - GHOST);
      vy_prim_pv_(i, size_y_ - j) = 0.0;
      pressure_prim_pv_(i, size_y_ - j)
                       = pressure_prim_pv_(i, size_y_ - GHOST);
    }
  }

  // sia_x
  for(auto j = 0; j < GHOST; ++j){
    for(auto i = 0; i < size_x_ + 1; ++i){
      density_prim_sia_x_(i, j) = density_prim_sia_x_(i, GHOST);
      vx_prim_sia_x_(i, j) = vx_prim_sia_x_(i, GHOST);
      vy_prim_sia_x_(i, j) = 0.0;
      pressure_prim_sia_x_(i, j) = pressure_prim_sia_x_(i, GHOST);

      density_prim_sia_x_(i, size_y_ - j - 1) =
          density_prim_sia_x_(i, size_y_ - GHOST - 1);
      vx_prim_sia_x_(i, size_y_ - j - 1) =
          vx_prim_sia_x_(i, size_y_ - GHOST - 1);
      vy_prim_sia_x_(i, size_y_ - j - 1) = 0.0;
      pressure_prim_sia_x_(i, size_y_ - j - 1) =
          pressure_prim_sia_x_(i, size_y_ - GHOST - 1);
    }
  }

  // sia_y
  for(auto j = 0; j < GHOST + 1; ++j){
    for(auto i = 0; i < size_x_; ++i){
      density_prim_sia_y_(i, j) = density_prim_sia_y_(i, GHOST);
      vx_prim_sia_y_(i, j) = vx_prim_sia_y_(i, GHOST);
      vy_prim_sia_y_(i, j) = 0.0;
      pressure_prim_sia_y_(i, j) = pressure_prim_sia_y_(i, GHOST);

      density_prim_sia_y_(i, size_y_ - j) =
          density_prim_sia_y_(i, size_y_ - GHOST);
      vx_prim_sia_y_(i, size_y_ - j) =
          vx_prim_sia_y_(i, size_y_ - GHOST);
      vy_prim_sia_y_(i, size_y_ - j) = 0.0;
      pressure_prim_sia_y_(i, size_y_ - j) =
          pressure_prim_sia_y_(i, size_y_ - GHOST);
    }
  }

  // via
  for(auto j = 0; j < GHOST; ++j){
    for(auto i = 0; i < size_x_; ++i){
      density_prim_via_(i, j) = density_prim_via_(i, GHOST);
      vx_prim_via_(i, j) = vx_prim_via_(i, GHOST);
      vy_prim_via_(i, j) = 0.0;
      pressure_prim_via_(i, j) = pressure_prim_via_(i, GHOST);

      density_prim_via_(i, size_y_ - j - 1) =
          density_prim_via_(i, size_y_ - GHOST - 1);
      vx_prim_via_(i, size_y_ - j - 1) =
          vx_prim_via_(i, size_y_ - GHOST - 1);
      vy_prim_via_(i, size_y_ - j - 1) = 0.0;
      pressure_prim_via_(i, size_y_ - j - 1) =
          pressure_prim_via_(i, size_y_ - GHOST - 1);
    }
  }
}

void Fluid2D::UpdateValue() noexcept
{
  // PV
  swap(density_prim_pv_, density_prim_pv_new_);
  swap(vx_prim_pv_, vx_prim_pv_new_);
  swap(vy_prim_pv_, vy_prim_pv_new_);
  swap(pressure_prim_pv_, pressure_prim_pv_new_);

  // SIA x
  swap(density_prim_sia_x_, density_prim_sia_x_new_);
  swap(vx_prim_sia_x_, vx_prim_sia_x_new_);
  swap(vy_prim_sia_x_, vy_prim_sia_x_new_);
  swap(pressure_prim_sia_x_, pressure_prim_sia_x_new_);

  // SIA y
  swap(density_prim_sia_y_, density_prim_sia_y_new_);
  swap(vx_prim_sia_y_, vx_prim_sia_y_new_);
  swap(vy_prim_sia_y_, vy_prim_sia_y_new_);
  swap(pressure_prim_sia_y_, pressure_prim_sia_y_new_);

  // VIA
  swap(density_prim_via_, density_prim_via_new_);
  swap(vx_prim_via_, vx_prim_via_new_);
  swap(vy_prim_via_, vy_prim_via_new_);
  swap(pressure_prim_via_, pressure_prim_via_new_);
}

void Fluid2D::Output(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << '-' << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto j = GHOST; j < size_y_ - GHOST; ++j){
    for(auto i = GHOST; i < size_x_ - GHOST; ++i){
      file << (x_[i] + x_[i + 1]) * 0.5 << ' '
           << (y_[j] + y_[j + 1]) * 0.5 << ' '
           << density_prim_via_(i, j) << ' '
           << vx_prim_via_(i, j) << ' '
           << vy_prim_via_(i, j) << ' '
           << pressure_prim_via_(i, j) << '\n';
    }
    file << std::endl;
  }

  file.close ();
}

void Fluid2D::OutputSiaX(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << "-SIA_X-" << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto j = GHOST; j < size_y_  - GHOST; ++j){
    for(auto i = GHOST; i < size_x_ + 1 - GHOST; ++i){
      file << x_[i] << ' '
           << (y_[j] + y_[j + 1]) * 0.5<< ' '
           << density_prim_sia_x_(i, j) << ' '
           << vx_prim_sia_x_(i, j) << ' '
           << vy_prim_sia_x_(i, j) << ' '
           << pressure_prim_sia_x_(i, j) << '\n';
    }
    file << std::endl;
  }

  file.close ();
}

void Fluid2D::OutputSiaY(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << "-SIA_Y-" << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto j = GHOST; j < size_y_  + 1 - GHOST; ++j){
    for(auto i = GHOST; i < size_x_ - GHOST; ++i){
      file << (x_[i] + x_[i + 1]) * 0.5 << ' '
           << y_[j] << ' '
           << density_prim_sia_y_(i, j) << ' '
           << vx_prim_sia_y_(i, j) << ' '
           << vy_prim_sia_y_(i, j) << ' '
           << pressure_prim_sia_y_(i, j) << '\n';
    }
    file << std::endl;
  }

  file.close ();
}
