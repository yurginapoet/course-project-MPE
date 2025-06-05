#include "header.h"

double f(nd &n, double t, int flag)
{
  switch (flag)
  {
  case 1:
    return 1;
  case 2:
    return 2 * t;
  case 3:
    return 3 * t * t;
  case 4:
    return 4 * t * t * t;
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
    return t;
    // return n.r * n.r + n.z + t;
  case 2:
    return t * t;
  case 3:
    return t * t * t;
  case 4:
    return t * t * t * t;
  case 5:
    return cos(t);

  default:
    return 0;
    // to be continued...
  }
}

double uInit(nd &n, double t, int flag) { return u(n, t, flag); }
