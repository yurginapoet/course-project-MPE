#pragma once
#include "func.h"

class Mesh
{
public:
  Mesh()
      : nodes(), blocks(), first(0), second(0), b(), u(), fictnum(), area(),
        matrix(), num(0), nx(0), ny(0), m(0), hx(0.0), hy(0.0), kx(0.0),
        ky(0.0), lambda(0.0), gamma(0.0), eps(0.0), w(0.0), iternum(0)
  {
  }

  void Clear();

  vector<vector<d>> blocks;

  vector<Point> nodes;
  vector<d> b, u;
  vector<int> first, second;
  vector<int> fictnum;
  Rect area;
  vector<vector<d>> matrix;
  int num, nx, ny, m, N;
  d hx, hy, kx, ky, rx, ry;
  d lambda, gamma;
  d eps, w;
  int iternum, testnum;

  void Out();
  void Input(string fin, string flg, string fparam);
  void CalcStep();
  void Nodes();
  void FindFictive();
  void Boundary(string ffirst, string fsecond);
  void RemoveSecondConditionals();
  void FillMatrix();

  bool IsFictive(int k);
  bool IsSecond(int k);
  bool IsFirst(int k);
  bool IsInnerNode(int i, int j, int k);

  void HaussZeidel(string filename);
  d Iteration(int i);
  void Results(string filename);
};