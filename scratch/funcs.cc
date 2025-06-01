#include "header.h"

double f(nd &n, double t, int flag)
{
  switch (flag)
  {
  case 1:
    return n.r * n.z + t;
  default:
    return 0;
    // to be continued...
  }
}

double u(nd &n, double t, int flag)
{
  switch (flag)
  {
  case 1:
    return n.r * n.z + t;
  default:
    return 0;
    // to be continued...
  }
}

double uInit(nd &n, double t, int flag)
{
  switch (flag)
  {
  case 1:
    return n.r * n.z + t;
  default:
    return 0;
    // to be continued...
  }
}
