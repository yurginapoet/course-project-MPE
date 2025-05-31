#include "header.h"

// Main function
int main()
{
  vector<nd> mesh;
  vector<el> elList;
  double gamma;
  int fnum;

  SLAE slae{}, LU{};
  LOS v{};
  TimeMesh timemesh{};

  input(mesh, elList, gamma, fnum);
  readTimeMesh(timemesh);
  readSplitTimeMesh(timemesh);

  GetPortraitSparseMatrix(mesh, elList, slae);

  return 0;
}