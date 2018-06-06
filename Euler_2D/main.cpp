#include <string>

#include "global_definitions.h"
#include "Fluid2D.h"
#include "FluidInitialize.h"

using namespace std;

int main()
{
  string fname_base("test");
  FluidInitialize init_func;

  Fluid2D fluid2d(64, 64, 0.0, 1.0, 0.0, 1.0, 1.4, init_func);

  fluid2d.MainLoop(2000, 4, fname_base);

  return 0;
}
