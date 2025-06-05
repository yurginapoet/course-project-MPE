#include "header.h"

double f(nd &n, double t, int flag)
{
  switch (flag)
  {
  case 1:
    // return 1;
    return 1 - 0;
  case 2:
    // return 2 * t;
    return 2 * t - 0;
  case 3:
    // return 3 * t * t;
    return 3 * t * t - 0;
  case 4:
    // return 4 * t * t * t;
    return 4 * t * t * t - 0;
  case 5:
    return -sin(t);
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
    // return t;
    return n.z + t;
  case 2:
    // return t * t;
    return n.z + t * t;
  case 3:
    // return t * t * t;
    return n.z + t * t * t;
  case 4:
    // return t * t * t * t;
    return n.z + t * t * t * t;
  case 5:
    return cos(t);

  default:
    return 0;
    // to be continued...
  }
}
