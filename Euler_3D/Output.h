template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
::Output(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << '-' << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto k = GHOST; k < nz_ - GHOST; ++k){
    for(auto j = GHOST; j < ny_ - GHOST; ++j){
      for(auto i = GHOST; i < nx_ - GHOST; ++i){
        file << (x_(i) + x_(i + 1)) * 0.5 << ' '
             << (y_(j) + y_(j + 1)) * 0.5 << ' '
             << (z_(k) + z_(k + 1)) * 0.5 << ' '
             << via_.d(i, j, k) << ' '
             << via_.u(i, j, k) << ' '
             << via_.v(i, j, k) << ' '
             << via_.w(i, j, k) << ' '
             << via_.p(i, j, k) << '\n';
      }
      file << std::endl;
    }
    file << std::endl;
  }

  file.close ();
}

template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
::OutputSliceX(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << "-SliceX-" << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  auto i = nx_ / 2;
  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto k = GHOST; k < nz_ - GHOST; ++k){
    for(auto j = GHOST; j < ny_ - GHOST; ++j){
      file << (y_(j) + y_(j + 1)) * 0.5 << ' '
           << (z_(k) + z_(k + 1)) * 0.5 << ' '
           << via_.d(i, j, k) << ' '
           << via_.u(i, j, k) << ' '
           << via_.v(i, j, k) << ' '
           << via_.w(i, j, k) << ' '
           << via_.p(i, j, k) << '\n';
    }
    file << std::endl;
  }

  file.close ();
}

template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
::OutputSliceY(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << "-SliceY-" << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  auto j = ny_ / 2;
  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto k = GHOST; k < nz_ - GHOST; ++k){
    for(auto i = GHOST; i < nx_ - GHOST; ++i){
      file << (x_(i) + x_(i + 1)) * 0.5 << ' '
           << (z_(k) + z_(k + 1)) * 0.5 << ' '
           << via_.d(i, j, k) << ' '
           << via_.u(i, j, k) << ' '
           << via_.v(i, j, k) << ' '
           << via_.w(i, j, k) << ' '
           << via_.p(i, j, k) << '\n';
    }
    file << std::endl;
  }

  file.close ();
}

template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
::OutputSliceZ(const std::string &filename, Int32 step) const
{
  std::stringstream out_count;
  out_count << "-SliceZ-" << std::setfill('0') << std::setw(4) << step;
  out_count << ".dat";

  auto k = nz_ / 2;
  std::ofstream file(filename + out_count.str(), std::ios::out);
  file << "# time = " << time_ << std::endl;
  for(auto j = GHOST; j < ny_ - GHOST; ++j){
    for(auto i = GHOST; i < nx_ - GHOST; ++i){
      file << (x_(i) + x_(i + 1)) * 0.5 << ' '
           << (y_(j) + y_(j + 1)) * 0.5 << ' '
           << via_.d(i, j, k) << ' '
           << via_.u(i, j, k) << ' '
           << via_.v(i, j, k) << ' '
           << via_.w(i, j, k) << ' '
           << via_.p(i, j, k) << '\n';
    }
    file << std::endl;
  }

  file.close ();
}
