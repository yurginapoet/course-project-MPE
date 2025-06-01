#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdio.h>
#include <vector>

typedef double d;
using namespace std;

struct Point
{
  d x;
  d y;
};

struct Rect
{
  d xmin, xmax, ymin, ymax;
};

d Step(int n, d k, d max, d min);

d F(d x, d y, int flag);
d U(d x, d y, int flag);
d Theta(d x, d y, int flag);
d Norm(int n, const vector<d> &v);
