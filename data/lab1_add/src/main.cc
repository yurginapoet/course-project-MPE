#include "Mesh.h"
#include "func.h"
#include <string>

int main()
{
  Mesh mesh;

  string fin = "data/in.txt";
  string ffirst = "data/first.txt";
  string fsecond = "data/second.txt";
  string flg = "data/lg.txt";
  string fparam = "data/param.txt";
  mesh.testnum = 1;

  cout << "0" << endl;
  mesh.Input(fin, flg, fparam);
  cout << "1" << endl;
  mesh.CalcStep();
  cout << 2 << endl;
  mesh.Nodes();
  cout << 3 << endl;
  mesh.FindFictive();
  cout << 4 << endl;
  mesh.Boundary(ffirst, fsecond);
  cout << 5 << endl;
  mesh.RemoveSecondConditionals();
  cout << 6 << endl;
  mesh.Out();
  mesh.FillMatrix();
  cout << 7 << endl;
  mesh.Out();
  mesh.HaussZeidel("data/results/hzz.txt");
  mesh.Results("data/results/test_approx_unregular_2nd_.txt");
  mesh.Clear();

  return 0;
}