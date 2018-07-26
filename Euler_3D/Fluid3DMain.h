template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
            ::MainLoop(Int32 loop_max, Int32 out_step,
                       const std::string &fname_base)
{
  Int32 output = 0;

  BoundaryConditions();

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

template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>::Update()
{
  Advection();
  BoundaryConditions();
  time_ += dt_;
  NewDt();
}

template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>::UpdateValue()
{
  // PV
  swap(pv_, pv_new_);

  // LIA x
  swap(lia_x_, lia_x_new_);

  // LIA y
  swap(lia_y_, lia_y_new_);

  // LIA z
  swap(lia_z_, lia_z_new_);

  // SIA x
  swap(sia_x_, sia_x_new_);

  // SIA y
  swap(sia_y_, sia_y_new_);

  // SIA z
  swap(sia_z_, sia_z_new_);

  // VIA
  swap(via_, via_new_);
}
