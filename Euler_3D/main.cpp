#include <string>

#include "global_definitions.h"
#include "Fluid3D.h"
#include "FluidInitialize.h"
#include "OutflowBoundary.h"
#include "ReflectBoundary.h"

using namespace std;

int main()
{
  string fname_base("test");
  FluidInitialize init_func;

  //              nx, ny, nz, xmin,xmax,ymin,ymax,zmin,zmax,gamma,init_func
  Fluid3D fluid3d(32, 32, 32, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.4, init_func);

  fluid3d.SetBCXLeft(OutflowBoundaryXLeft);
  fluid3d.SetBCXRight(OutflowBoundaryXRight);

  fluid3d.SetBCYLeft(OutflowBoundaryYLeft);
  fluid3d.SetBCYRight(OutflowBoundaryYRight);

  fluid3d.SetBCZLeft(OutflowBoundaryZLeft);
  fluid3d.SetBCZRight(OutflowBoundaryZRight);

  fluid3d.MainLoop(1600, 4, fname_base);

  return 0;
}
