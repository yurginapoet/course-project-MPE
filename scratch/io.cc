#include "header.h"

void input(vector<nd> &mesh, vector<el> &elList, double &gamma, int &fnum)
{
  int n, k, node;
  double r, z, l;

  // Read nodes
  ifstream in("../data/nds.txt");
  in >> n;
  mesh.resize(n);
  for (int i = 0; i < n; i++)
  {
    in >> r >> z;
    mesh[i] = nd(i, r, z);
  }
  in.close();

  // Read elements
  in.open("../data/els.txt");
  in >> n;
  elList.resize(n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      in >> node;
      elList[i].nds[j] = mesh[node];
    }
  }
  in.close();

  // Read lambda
  in.open("../data/lambda.txt");
  for (int i = 0; i < n; i++)
    in >> elList[i].lambda;
  in.close();

  // Read sigma
  in.open("../data/gamma.txt");
  in >> gamma;
  in.close();

  // Read fnum
  in.open("../data/fnum.txt");
  in >> fnum;
  in.close();
}

void readTimeMesh(TimeMesh &time)
{
  int countTimeLayer = 0;

  ifstream timeMesh("../data/timeMesh.txt");

  timeMesh >> countTimeLayer;
  time.t.resize(countTimeLayer);

  for (int ti = 0; ti < countTimeLayer; ti++)
    timeMesh >> time.t[ti];

  timeMesh.close();
}

void readSplitTimeMesh(TimeMesh &time)
{
  auto &tNew = time.tNew;

  const int tSize = time.t.size();

  ifstream splittingTimeMesh("../data/timeMeshSplit.txt");

  int nk = 0;
  tNew.resize(1, time.t[0]);

  for (int ti = 0, j = 1; ti < tSize - 1; ti++, j++)
  {
    int countIntervals = 0;

    double coef = 0;
    double step = 0;

    splittingTimeMesh >> countIntervals >> coef;
    nk += countIntervals;
    tNew.resize(nk + 1);

    if (coef != 1)
    {
      double sumProgression = (pow(coef, countIntervals) - 1.) / (coef - 1.);
      step = (time.t[ti + 1] - time.t[ti]) / sumProgression;

      int jk = 1;
      for (j; j < nk; j++, jk++)
        tNew[j] = time.t[ti] + step * (pow(coef, jk) - 1.) / (coef - 1.);
    }
    else
    {
      step = (time.t[ti + 1] - time.t[ti]) / countIntervals;

      int jk = 1;
      for (j; j < nk; j++, jk++)
        tNew[j] = time.t[ti] + step * jk;
    }

    tNew[j] = time.t[ti + 1];
  }

  time.q.resize(time.tNew.size());
}

void readBoundaryCondition(vector<BoundaryConditions> &cond)
{
  int numEdgeConditions = 0;

  ifstream conditions("../data/boundaryConditions.txt");
  conditions >> numEdgeConditions;
  cond.resize(numEdgeConditions);

  for (int i = 0; i < numEdgeConditions; i++)
  {
    auto &VerticesNumbers = cond[i].VerticesNumbers;
    int numVertex = 0, type = 0, function = 0;

    conditions >> type >> function >> numVertex;

    cond[i].type = type;
    cond[i].function = function - 1;
    VerticesNumbers.resize(numVertex);

    for (int k = 0; k < numVertex; k++)
      conditions >> VerticesNumbers[k];
  }
}