#include "header.h"

// Main function
int main()
{
  vector<nd> mesh;
  vector<el> elList;
  double sigma;
  int fnum;

  int flag = 1;

  SLAE slae{}, LU{};
  LOS v{};
  TimeMesh timemesh{};
  vector<BoundaryConditions> conds{};

  input(mesh, elList, sigma, fnum);
  readTimeMesh(timemesh);
  readSplitTimeMesh(timemesh);
  readBoundaryCondition(conds);

  GetPortraitSparseMatrix(mesh, elList, slae);

  Solve(mesh, elList, timemesh, slae, conds, flag, sigma);

  return 0;
}