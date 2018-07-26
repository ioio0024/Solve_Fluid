#include <string>

#include "global_definitions.h"
#include "CipCsl3.h"
#include "Fluid3D.h"
#include "FluidInitialize.h"
#include "OutflowBoundary.h"
#include "ReflectBoundary.h"

using namespace std;

int main()
{
  string fname_base("test");
  FluidInitialize init_func;

  Fluid3D<OutflowBoundaryXLeft, OutflowBoundaryXRight,
          OutflowBoundaryYLeft, OutflowBoundaryYRight,
          OutflowBoundaryZLeft, OutflowBoundaryZRight,
          CipCsl3X, CipCsl3Y, CipCsl3Z>
          fluid3d(32, 32, 32, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.4, init_func);
  //              nx, ny, nz, xmin,xmax,ymin,ymax,zmin,zmax,gamma,init_func

  fluid3d.MainLoop(1600, 4, fname_base);

  return 0;
}
