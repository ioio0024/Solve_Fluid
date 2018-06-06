#include "global_definitions.h"
#include "Fluid3D.h"

#include <fstream>
#include <iomanip>
#include <sstream>

void Fluid3D::Output(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << '-' << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto k = GHOST; k < nz_ - GHOST; ++k){
    for(auto j = GHOST; j < ny_ - GHOST; ++j){
      for(auto i = GHOST; i < nx_ - GHOST; ++i){
        file << (x_[i] + x_[i + 1]) * 0.5 << ' '
             << (y_[j] + y_[j + 1]) * 0.5 << ' '
             << (z_[k] + z_[k + 1]) * 0.5 << ' '
             << density_via_(i, j, k) << ' '
             << vx_via_(i, j, k) << ' '
             << vy_via_(i, j, k) << ' '
             << vz_via_(i, j, k) << ' '
             << pressure_via_(i, j, k) << '\n';
      }
      file << std::endl;
    }
    file << std::endl;
  }

  file.close ();
}

void Fluid3D::OutputSliceX(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << "-SliceX-" << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  auto i = nx_ / 2;
  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto k = GHOST; k < nz_ - GHOST; ++k){
    for(auto j = GHOST; j < ny_ - GHOST; ++j){
      file << (y_[j] + y_[j + 1]) * 0.5 << ' '
           << (z_[k] + z_[k + 1]) * 0.5 << ' '
           << density_via_(i, j, k) << ' '
           << vx_via_(i, j, k) << ' '
           << vy_via_(i, j, k) << ' '
           << vz_via_(i, j, k) << ' '
           << pressure_via_(i, j, k) << '\n';
    }
    file << std::endl;
  }

  file.close ();
}

void Fluid3D::OutputSliceY(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << "-SliceY-" << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  auto j = ny_ / 2;
  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto k = GHOST; k < nz_ - GHOST; ++k){
    for(auto i = GHOST; i < nx_ - GHOST; ++i){
      file << (x_[i] + x_[i + 1]) * 0.5 << ' '
           << (z_[k] + z_[k + 1]) * 0.5 << ' '
           << density_via_(i, j, k) << ' '
           << vx_via_(i, j, k) << ' '
           << vy_via_(i, j, k) << ' '
           << vz_via_(i, j, k) << ' '
           << pressure_via_(i, j, k) << '\n';
    }
    file << std::endl;
  }

  file.close ();
}

void Fluid3D::OutputSliceZ(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << "-SliceZ-" << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  auto k = nz_ / 2;
  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto j = GHOST; j < ny_ - GHOST; ++j){
    for(auto i = GHOST; i < nx_ - GHOST; ++i){
      file << (x_[i] + x_[i + 1]) * 0.5 << ' '
           << (y_[j] + y_[j + 1]) * 0.5 << ' '
           << density_via_(i, j, k) << ' '
           << vx_via_(i, j, k) << ' '
           << vy_via_(i, j, k) << ' '
           << vz_via_(i, j, k) << ' '
           << pressure_via_(i, j, k) << '\n';
    }
    file << std::endl;
  }

  file.close ();
}
