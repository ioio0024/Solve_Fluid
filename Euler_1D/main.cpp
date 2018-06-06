#include <string>

#include "global_definitions.h"
#include "Fluid1D.h"
#include "FluidInitialize.h"

using namespace std;

int main()
{
  string fname_base("test");

  Fluid1D fluid1d(200, 0.0, 1.0, 1.4, FluidInitialize);

  fluid1d.MainLoop(500, 1, fname_base);

  return 0;
}
