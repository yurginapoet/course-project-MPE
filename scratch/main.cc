// #include "header.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>
using namespace std;

// структура узла сетки
struct nd
{
  int gl_num;
  double r, z;
  nd(int _num, double _r, double _z)
  {
    gl_num = _num;
    r = _r;
    z = _z;
  }
  nd() = default;
};

// структура элемента
struct el
{
  double lambda;
  nd nds[3];
};
// structure for time layer?????????????
struct timelayer
{
  double t;      // value of current timestep
  int intervals; // number of intervals on current timestep
  double coef;   // coefficient of intervals on current time step
};

// time mesh
struct TimeMesh
{
  vector<double> t{};
  vector<double> tNew{};

  vector<vector<double>> q{};
  vector<double> qti{};
  vector<double> qti_1{};
  vector<double> qti_2{};

  void swap();
  void saveWeights(int ti);
};

// time mesh funcs
void TimeMesh::swap()
{
  std::swap(qti_2, qti_1);
  std::swap(qti_1, qti);
}
void TimeMesh::saveWeights(int ti) { q[ti] = qti; }

// sparse matrix structure
struct SparseMatrix
{
  vector<int> ig{}, jg{};
  vector<double> ggl{}, ggu{}, di{};
};

// SLAE structure
struct SLAE
{
  SparseMatrix A{};
  vector<double> q{}, b{};
};

// LOS structure
struct LOS
{
  vector<double> r1{}, rk{}, z1{}, p1{}, Ar{}, p{}, mult{};
};

struct BoundaryConditions
{
  int type{}, function{};
  vector<int> VerticesNumbers{};
};

void GetPortraitSparseMatrix(vector<nd> &mesh, vector<el> &elList, SLAE &slae)
{
  auto &nodeCoord = mesh;
  auto &elements = elList;

  auto &ig = slae.A.ig, &jg = slae.A.jg;
  auto &b = slae.b;

  vector<std::set<int>> rowCount{};

  const int sizeSlae = nodeCoord.size();

  slae.A.di.resize(sizeSlae);
  slae.b.resize(sizeSlae);
  slae.q.resize(sizeSlae);
  ig.resize(sizeSlae + 1);
  rowCount.resize(sizeSlae);

  for (auto &elem : elements)
  {
    for (auto &i : elem.nds)
    {
      for (auto &j : elem.nds)
        if (j.gl_num < i.gl_num)
          rowCount[i.gl_num].insert(j.gl_num);
    }
  }

  ig[0] = 0;
  for (int i = 0; i < sizeSlae; i++)
  {
    ig[i + 1] = ig[i] + rowCount[i].size();
    for (auto j : rowCount[i])
      jg.push_back(j);
  }

  slae.A.ggl.resize(ig[sizeSlae]);
  slae.A.ggu.resize(ig[sizeSlae]);
}

void localOptimalSchemeLU(SLAE &slae, SLAE &LU, LOS &v, int maxIter, double eps)
{
  auto &r1 = v.r1, &z1 = v.z1, &p1 = v.p1, &mult = v.mult, &rk = v.rk,
       &Ar = v.Ar, &p = v.p;
  auto &q = slae.q;
  const int n = slae.A.di.size();
  double normb = 0;

  p.resize(n);
  r1.resize(n);
  z1.resize(n);
  p1.resize(n);
  mult.resize(n);
  rk.resize(n);
  Ar.resize(n);

  calcDiscrepancy(slae, v, q, normb);
  calcY(LU, r1, r1);
  calcX(LU, r1, z1);
  multOfMatrix(slae.A, z1, p1);
  calcY(LU, p1, p1);
  double scalarr = scalarMult(r1, r1), discrepancy = sqrt(scalarr / normb);

  for (int k = 1; k < maxIter && discrepancy > eps; k++)
  {
    double scalarp = scalarMult(p1, p1), alpha = scalarMult(p1, r1) / scalarp;
    calcVectorMultCoef(z1, alpha, mult);
    calcSumVectors(q, v.mult, q);
    calcVectorMultCoef(p1, -alpha, mult);
    calcSumVectors(r1, mult, r1);

    calcX(LU, r1, rk);
    multOfMatrix(slae.A, rk, Ar);
    calcY(LU, Ar, p);
    double betta = -scalarMult(p1, p) / scalarp;
    calcVectorMultCoef(z1, betta, mult);
    calcSumVectors(rk, mult, z1);
    calcVectorMultCoef(p1, betta, mult);
    calcSumVectors(p, mult, p1);
    discrepancy = sqrt(scalarMult(r1, r1) / scalarr);
    std::cout << k << " " << discrepancy << std::endl;
  }
  normb = 0;
  calcDiscrepancy(slae, v, q, normb);
  discrepancy = sqrt(scalarMult(r1, r1) / normb);
  std::cout << "Final discrepancy: " << discrepancy << std::endl;
}

void calcLU(SLAE &slae, SLAE &LU)
{
  auto &ig = slae.A.ig, &jg = slae.A.jg;
  auto &ggl = slae.A.ggl, &ggu = slae.A.ggu, &di = slae.A.di, &L = LU.A.ggl,
       &U = LU.A.ggu, &diL = LU.A.di;
  LU.b = slae.b;
  LU.A.ig = ig;
  LU.A.jg = jg;

  const int sizeSlae = di.size(), sizeTriangles = ig[sizeSlae];

  diL.resize(sizeSlae);
  L.resize(sizeTriangles);
  U.resize(sizeTriangles);

  for (int i = 0; i < sizeSlae; i++)
  {
    double sumDi = 0;
    int i0 = ig[i];
    int i1 = ig[i + 1];

    for (int k = i0; k < i1; k++)
    {
      double suml = 0, sumu = 0;
      int j = jg[k];
      int j0 = ig[j];
      int j1 = ig[j + 1];

      for (int ik = i0, kj = j0; ik < i1 && kj < j1;)
      {
        if (jg[ik] > jg[kj])
          kj++;
        else if (jg[ik] < jg[kj])
          ik++;
        else
        {
          suml += L[ik] * U[kj];
          sumu += L[kj] * U[ik];
          ik++;
          kj++;
        }
      }

      L[k] = (ggl[k] - suml);
      U[k] = (ggu[k] - sumu) / diL[j];
      sumDi += L[k] * U[k];
    }
    diL[i] = di[i] - sumDi;
  }
}

void calcY(SLAE &LU, vector<double> &b, vector<double> &y)
{
  auto &ig = LU.A.ig, &jg = LU.A.jg;
  auto &di = LU.A.di, &L = LU.A.ggl;
  const int sizeSlae = di.size();

  for (int i = 0; i < sizeSlae; i++)
  {
    double sum = 0;
    int i0 = ig[i], i1 = ig[i + 1];

    for (i0; i0 < i1; i0++)
    {
      int j = jg[i0];
      sum += L[i0] * y[j];
    }

    y[i] = (b[i] - sum) / di[i];
  }
}

void calcX(SLAE &LU, vector<double> &y, vector<double> &x)
{
  auto &ig = LU.A.ig, &jg = LU.A.jg;
  auto &U = LU.A.ggu;
  const int sizeSlae = LU.A.di.size();
  vector<double> v = y;

  for (int i = sizeSlae - 1; i >= 0; i--)
  {
    x[i] = v[i];
    int i0 = ig[i], i1 = ig[i + 1];

    for (i0; i0 < i1; i0++)
    {
      int j = jg[i0];
      v[j] -= x[i] * U[i0];
    }
  }
}

void multOfMatrix(SparseMatrix &A, vector<double> &x, vector<double> &F)
{
  auto &ig = A.ig, &jg = A.jg;
  auto &di = A.di, &ggl = A.ggl, &ggu = A.ggu;
  const int sizeA = di.size();

  for (int i = 0; i < sizeA; i++)
  {
    F[i] = di[i] * x[i];
    int i0 = ig[i], i1 = ig[i + 1];

    for (i0; i0 < i1; i0++)
    {
      int j = jg[i0];
      F[i] += ggl[i0] * x[j];
      F[j] += ggu[i0] * x[i];
    }
  }
}

void calcDiscrepancy(SLAE &slae, LOS &v, vector<double> &x, double &normb)
{
  auto &ig = slae.A.ig, &jg = slae.A.jg;
  auto &ggl = slae.A.ggl, &ggu = slae.A.ggu, &di = slae.A.di, &b = slae.b,
       &r1 = v.r1;
  const int sizeSlae = di.size();

  for (int i = 0; i < sizeSlae; i++)
  {
    normb += b[i] * b[i];
    r1[i] = b[i] - di[i] * x[i];
    int i0 = ig[i], i1 = ig[i + 1];
    for (i0; i0 < i1; i0++)
    {
      int j = jg[i0];
      r1[i] -= ggl[i0] * x[j];
      r1[j] -= ggu[i0] * x[i];
    }
  }
}

void calcVectorMultCoef(vector<double> &a, double coef, vector<double> &res)
{
  const int n = a.size();

  for (int i = 0; i < n; i++)
    res[i] = a[i] * coef;
}

void calcSumVectors(vector<double> &a, vector<double> &b, vector<double> &res)
{
  const int n = a.size();

  for (int i = 0; i < n; i++)
    res[i] = a[i] + b[i];
}

double scalarMult(vector<double> &a, vector<double> &b)
{
  int n = a.size();
  double res = 0;
  for (int i = 0; i < n; i++)
    res += a[i] * b[i];
  return res;
}

void clearSLAE(SLAE &slae)
{
  auto &A = slae.A;
  auto &b = slae.b;
  auto &q = slae.q;

  const int sizeSLAE = A.di.size();
  const int countNonZeroElems = A.ig[sizeSLAE];

  for (int i = 0; i < sizeSLAE; i++)
  {
    A.di[i] = 0;
    q[i] = 0;
    b[i] = 0;
  }

  for (int i = 0; i < countNonZeroElems; i++)
  {
    A.ggl[i] = 0;
    A.ggu[i] = 0;
  }
}

void input(vector<nd> &mesh, vector<el> &elList, double &gamma, int &fnum,
           vector<timelayer> &time)
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

  // Read gamma
  in.open("../data/gamma.txt");
  in >> gamma;
  in.close();

  // Read fnum
  in.open("../data/fnum.txt");
  in >> fnum;
  in.close();

  // Read time steps values
  in.open("../data/time.txt");
  in >> k;
  time.resize(k);
  for (int i = 0; i < k; i++)
    in >> time[i].t;
  in.close();

  // Read time number of intervals & coefficients
  in.open("../data/timeCoef.txt");
  for (int i = 0; i < k; i++)
    in >> time[i].intervals >> time[i].coef;
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

// solver
void Solve(vector<nd> &mesh, vector<el> &elList, TimeMesh &timemesh, SLAE &slae,
           vector<BoundaryConditions> &conds)
{
  const int tSize = timemesh.tNew.size();

  getWeightsInitU(timemesh, mesh);
  timemesh.q[0] = timemesh.qti_2;
  timemesh.q[1] = timemesh.qti_1;

  for (int ti = 2; ti < tSize; ti++)
  {
    SLAE LU{};
    LOS vectors{};

    GetGlobalMatrixAndVector(mesh, time, slae, conds, ti);

    calcLU(slae, LU);
    localOptimalSchemeLU(slae, LU, vectors, 10000, 1e-14);

    timemesh.qti = slae.q;
    timemesh.saveWeights(ti);
    timemesh.swap();

    clearSLAE(slae);
  }
}
// init q0 q1
void getWeightsInitU(TimeMesh &timemesh, vector<nd> &mesh)
{
  auto &qti = timemesh.qti;
  auto &qti_1 = timemesh.qti_1;
  auto &qti_2 = timemesh.qti_2;

  const int sizeMatrix = mesh.size();

  qti.resize(sizeMatrix);
  qti_1.resize(sizeMatrix);
  qti_2.resize(sizeMatrix);

  for (int i = 0; i < sizeMatrix; i++)
  {
    qti_2[i] = uInit(mesh[i], timemesh.tNew[0]);
    qti_1[i] = uInit(mesh[i], timemesh.tNew[1]);
  }
}
// add switch!!!!!!!!1
double uInit(nd &n, double t) { return n.r * n.r + n.z * n.z + t * t; }

// boundaries input
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

void GetGlobalMatrixAndVector(vector<nd> &mesh, vector<el> elList,
                              TimeMesh &timemesh, SLAE &slae,
                              vector<BoundaryConditions> &cond, int ti)
{
  auto &elements = elList;
  auto &vertexCoord = mesh;

  auto &A = slae.A;
  auto &b = slae.b;

  double t = timemesh.tNew[ti];
  double deltaT = timemesh.tNew[ti] - timemesh.tNew[ti - 2];
  double deltaT0 = timemesh.tNew[ti] - timemesh.tNew[ti - 1];
  double deltaT1 = timemesh.tNew[ti - 1] - timemesh.tNew[ti - 2];

  const int sizeSlae = vertexCoord.size();

  vector<vector<double>> M(3), G(3), tempMatrix(3);
  vector<double> locb(3, 0), tempVector(3, 0);

  // изменяем размер векторов
  for (int i = 0; i < 3; i++)
  {
    M[i].resize(3);
    G[i].resize(3);
    tempMatrix[i].resize(3);
  }

  for (auto &elem : elements)
  {
    // НУЖНО ДОБАВИТЬ СИГМУ
    getSigmaWeights(mesh.vertexCoord, elem, t);

    // ПРОСТО ТУПО ВЫЧИСЛЯЕМ ЛОКАЛЬНЫЕ МАТРИЦЫ И ВЕКТОРЫ
    getG(elem, G, t);
    getM(elem, M);
    getLocalb(vertexCoord, elem, locb, t);

    multiplyMatrixToCoef(M, (deltaT + deltaT0) / (deltaT * deltaT0),
                         tempMatrix);

    addLocalMatrixToGlobal(A, elem.localVertex, tempMatrix);
    addLocalMatrixToGlobal(A, elem.localVertex, G);

    addLocalVectorToGlobal(b, elem.localVertex, locb);

    multiplyMatrixToVector(M, timemesh.qti_2, tempVector, elem.localVertex);
    multiplyVectorToCoef(tempVector, -deltaT0 / (deltaT * deltaT1));
    addLocalVectorToGlobal(b, elem.localVertex, tempVector);

    multiplyMatrixToVector(M, timemesh.qti_1, tempVector, elem.localVertex);
    multiplyVectorToCoef(tempVector, deltaT / (deltaT1 * deltaT0));
    addLocalVectorToGlobal(b, elem.localVertex, tempVector);
  }

  addSecondBoundaryCondition(slae, cond, vertexCoord, t);
  addFirstBoundaryCondition(slae, cond, vertexCoord, t);
}

// Main function
int main()
{
  vector<nd> mesh;
  vector<el> elList;
  vector<timelayer> time; // рудимент видимо нужно будет убрать
  double gamma;
  int fnum;

  SLAE slae{}, LU{};
  LOS v{};
  TimeMesh timemesh{};
  vector<BoundaryConditions> conds{};

  input(mesh, elList, gamma, fnum, time);
  readTimeMesh(timemesh);
  readSplitTimeMesh(timemesh);

  GetPortraitSparseMatrix(mesh, elList, slae);

  return 0;
}