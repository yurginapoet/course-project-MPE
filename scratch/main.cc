#include "header.h"

// Main function
int main()
{
  vector<nd> mesh;
  vector<el> elList;
  double sigma;
  int flag;

  SLAE slae{}, LU{};
  LOS v{};
  TimeMesh timemesh{};
  vector<BoundaryConditions> conds{};

  input(mesh, elList, sigma, flag);
  input_timemesh(timemesh);
  input_split_timemesh(timemesh);
  input_boundary(conds);

  portrait(mesh, elList, slae);
  Solver(mesh, elList, timemesh, slae, conds, flag, sigma);
  print_solution(timemesh, mesh);

  return 0;
}