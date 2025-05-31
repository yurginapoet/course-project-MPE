#include "header.h"

// Main function
int main()
{
  vector<nd> mesh;
  vector<el> elList;
  double sigma;
  int fnum;

  SLAE slae{}, LU{};
  LOS v{};
  TimeMesh timemesh{};

  input(mesh, elList, sigma, fnum);
  readTimeMesh(timemesh);
  readSplitTimeMesh(timemesh);

  GetPortraitSparseMatrix(mesh, elList, slae);

  return 0;
}