#include "func.h"

d Step(int n, d k, d max, d min)
{
  d s = 0;
  for (int i = 0; i < n; i++)
    s += pow(k, i);
  return (max - min) / s;
}

d F(d x, d y, int flag)
{
  switch (flag)
  {
  case 1:
    return x + y;
  case 2:
    return y * y + x * x - 4;
  case 3:
    return pow(x, 3) + pow(y, 3) - 6 * (x + y);
  case 4:
    return pow(y, 4) + pow(x, 4) - 12 * (y * y + x * x);
  case 5:
    return pow(y, 4) - 12 * y * y;
  case 6:
    return cos(y);
  default:
    return 0;
  }
}

d U(d x, d y, int flag)
{
  switch (flag)
  {
  case 1:
    return x + y;
  case 2:
    return y * y + x * x;
  case 3:
    return pow(x, 3) + pow(y, 3);
  case 4:
    return pow(y, 4) + pow(x, 4);
  case 5:
    return pow(y, 4);
  case 6:
    return cos(y);
  default:
    return 0;
  }
}

d Theta(d x, d y, int flag)
{
  switch (flag)
  {
  case 1:
    return -1;
  case 2:
    return -2 * x;
  case 3:
    return -3 * x * x;
  case 4:
    return -4 * x * x * x;
  case 5:
    return 0;
  case 6:
    return 0;
  default:
    return 0;
  }
}

d Norm(int n, const vector<d> &v)
{
  d s = 0;
  for (int i = 0; i < n; i++)
    s += v[i] * v[i];
  return sqrt(s);
}

// // // degree: 4
// d F(d x, d y) { return pow(y, 4) - 12 * y * y; }
// d U(d x, d y) { return pow(y, 4); }

// // // degree: 2
// d F(d x, d y) { return y * y + x * x - 4; }
// d U(d x, d y) { return y * y + x * x; }
// d Theta(d x, d y) { return -2 * x; }

// // // degree: 1
// d F(d x, d y) { return x + y; }
// d U(d x, d y) { return x + y; }
// d Theta(d x, d y) { return -1; }

// // // degree: 3
// d F(d x, d y) { return pow(x, 3) + pow(y, 3) - 6 * (x + y); }
// d U(d x, d y) { return pow(x, 3) + pow(y, 3); }
// d Theta(d x, d y) { return -3 * x * x; }
