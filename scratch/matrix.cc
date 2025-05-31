#include "header.h"

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

/*=================================== LOS ====================================*/

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

/*======================= MATRIX & VECTOR OPERATIONS =========================*/
void addLocalMatrixToGlobal(SparseMatrix &A, vector<int> &localVertex,
                            vector<vector<double>> &localMatrix)
{
  const int sizeLocal = 6;

  for (int i = 0; i < sizeLocal; i++)
    for (int j = 0; j < sizeLocal; j++)
    {
      double elem = localMatrix[i][j];
      addElemToGlobalMatrix(A, localVertex[i], localVertex[j], elem);
    }
}

void addLocalVectorToGlobal(vector<double> &b, vector<int> &localVertex,
                            vector<double> &bLocal)
{
  const int sizeLocal = bLocal.size();

  for (int i = 0; i < sizeLocal; i++)
    b[localVertex[i]] += bLocal[i];
}

void addElemToGlobalMatrix(SparseMatrix &A, int i, int j, double elem)
{
  auto &ig = A.ig, &jg = A.jg;
  auto &ggl = A.ggl, &ggu = A.ggu, &di = A.di;

  if (i == j)
    di[i] += elem;
  else if (i > j)
  {
    int beg = ig[i], end = ig[i + 1] - 1;
    while (jg[beg] != j)
    {
      int ind = (beg + end) / 2;
      if (jg[ind] < j)
        beg = ind + 1;
      else
        end = ind;
    }
    ggl[beg] += elem;
  }
  else
  {
    int beg = ig[j], end = ig[j + 1] - 1;
    while (jg[beg] != i)
    {
      int ind = (beg + end) / 2;
      if (jg[ind] < i)
        beg = ind + 1;
      else
        end = ind;
    }
    ggu[beg] += elem;
  }
}

void multiplyMatrixToVector(vector<vector<double>> &matrix, vector<double> &vec,
                            vector<double> &result, vector<int> &localNum)
{
  const int sizeMatrix = matrix.size();

  for (int i = 0; i < sizeMatrix; i++)
  {
    double sum = 0;

    for (int j = 0; j < sizeMatrix; j++)
      sum += matrix[i][j] * vec[localNum[j]];

    result[i] = sum;
  }
}

void multiplyMatrixToCoef(vector<vector<double>> &matrix, double coef,
                          vector<vector<double>> &resultMatrix)
{
  const int sizeMatrix = matrix.size();

  for (int i = 0; i < sizeMatrix; i++)
    for (int j = 0; j < sizeMatrix; j++)
      resultMatrix[i][j] = coef * matrix[i][j];
}

void multiplyVectorToCoef(vector<double> &vector, double coef)
{
  const int size = vector.size();

  for (int i = 0; i < size; i++)
    vector[i] *= coef;
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

void calcVectorMultCoef(vector<double> &a, double coef, vector<double> &res)
{
  const int n = a.size();

  for (int i = 0; i < n; i++)
    res[i] = a[i] * coef;
}
