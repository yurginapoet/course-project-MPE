#include "header.h"

/*======================= MATRIX & VECTOR OPERATIONS =========================*/
void add_mx(SparseMatrix &A, el elem, vector<vector<double>> &localMatrix)
{
  const int sizeLocal = 3;
  auto locV = elem.nds;

  for (int i = 0; i < sizeLocal; i++)
    for (int j = 0; j < sizeLocal; j++)
    {
      double elem = localMatrix[i][j];
      add_el(A, locV[i].gl_num, locV[j].gl_num, elem);
    }
}

void add_vec(vector<double> &b, el &elem, vector<double> &bLocal)
{
  const int sizeLocal = bLocal.size();
  auto locV = elem.nds;

  for (int i = 0; i < sizeLocal; i++)
    b[locV[i].gl_num] += bLocal[i];
}

void add_el(SparseMatrix &A, int i, int j, double elem)
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

void mult_mx_vec(vector<vector<double>> &matrix, vector<double> &vec,
                 vector<double> &result, el &elem)
{
  const int sizeMatrix = matrix.size();
  auto locV = elem.nds;

  for (int i = 0; i < sizeMatrix; i++)
  {
    double sum = 0;

    for (int j = 0; j < sizeMatrix; j++)
      sum += matrix[i][j] * vec[locV[j].gl_num];

    result[i] = sum;
  }
}

void mult_mx_num(vector<vector<double>> &matrix, double coef,
                 vector<vector<double>> &resultMatrix)
{
  const int sizeMatrix = matrix.size();

  for (int i = 0; i < sizeMatrix; i++)
    for (int j = 0; j < sizeMatrix; j++)
      resultMatrix[i][j] = coef * matrix[i][j];
}

void mult_vec_num(vector<double> &vector, double coef)
{
  const int size = vector.size();

  for (int i = 0; i < size; i++)
    vector[i] *= coef;
}

void sum_vec(vector<double> &a, vector<double> &b, vector<double> &res)
{
  const int n = a.size();

  for (int i = 0; i < n; i++)
    res[i] = a[i] + b[i];
}

double dotproduct(vector<double> &a, vector<double> &b)
{
  int n = a.size();
  double res = 0;
  for (int i = 0; i < n; i++)
    res += a[i] * b[i];
  return res;
}

void mult_vec_num(vector<double> &a, double coef, vector<double> &res)
{
  const int n = a.size();

  for (int i = 0; i < n; i++)
    res[i] = a[i] * coef;
}

void mult_smx_vec(SparseMatrix &A, vector<double> &x, vector<double> &F)
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

void calc_discrepancy(SLAE &slae, LOS &v, vector<double> &x, double &normb)
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

void portrait(vector<nd> &mesh, vector<el> &elList, SLAE &slae)
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

void losLU(SLAE &slae, SLAE &LU, LOS &v, int maxIter, double eps)
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

  calc_discrepancy(slae, v, q, normb);
  if (std::abs(normb) < 1e-10)
  {
    std::cerr << "Error: normb is zero or too small: " << normb << std::endl;
    exit(1);
  }

  calcY(LU, r1, r1);
  calcX(LU, r1, z1);
  mult_smx_vec(slae.A, z1, p1);
  calcY(LU, p1, p1);
  double scalarr = dotproduct(r1, r1);
  double discrepancy = sqrt(scalarr / normb);
  if (std::isnan(discrepancy) || std::isinf(discrepancy))
  {
    std::cerr << "Error: Initial discrepancy is NaN or Inf" << std::endl;
    exit(1);
  }

  for (int k = 1; k < maxIter && discrepancy > eps; k++)
  {
    double scalarp = dotproduct(p1, p1);
    if (std::abs(scalarp) < 1e-35)
    {
      std::cerr << "Error: scalarp is zero or too small at iteration " << k
                << "it equals:" << scalarp << std::endl;
      exit(1);
    }
    double alpha = dotproduct(p1, r1) / scalarp;
    mult_vec_num(z1, alpha, mult);
    sum_vec(q, v.mult, q);
    mult_vec_num(p1, -alpha, mult);
    sum_vec(r1, mult, r1);

    calcX(LU, r1, rk);
    mult_smx_vec(slae.A, rk, Ar);
    calcY(LU, Ar, p);
    double betta = -dotproduct(p1, p) / scalarp;
    mult_vec_num(z1, betta, mult);
    sum_vec(rk, mult, z1);
    mult_vec_num(p1, betta, mult);
    sum_vec(p, mult, p1);
    discrepancy = sqrt(dotproduct(r1, r1) / normb);
    if (std::isnan(discrepancy) || std::isinf(discrepancy))
    {
      std::cerr << "Error: Discrepancy is NaN or Inf at iteration " << k
                << std::endl;
      exit(1);
    }
    std::cout << k << " " << discrepancy << std::endl;
  }
  normb = 0;
  calc_discrepancy(slae, v, q, normb);
  discrepancy = sqrt(dotproduct(r1, r1) / normb);
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

  diL.resize(sizeSlae, 0.0);
  L.resize(sizeTriangles, 0.0);
  U.resize(sizeTriangles, 0.0);

  for (int i = 0; i < sizeSlae; i++)
  {
    double sumDi = 0.0;
    int i0 = ig[i], i1 = ig[i + 1];

    // Вычисляем L[k] и U[k]
    for (int k = i0; k < i1; k++)
    {
      int j = jg[k];
      double suml = 0.0, sumu = 0.0;

      // Находим соответствующие элементы для суммирования
      for (int m = ig[j]; m < ig[j + 1]; m++)
      {
        for (int n = ig[i]; n < ig[i + 1]; n++)
        {
          if (jg[m] == jg[n] && jg[m] < j)
          {
            suml += L[n] * U[m];
            sumu += L[m] * U[n];
          }
        }
      }

      L[k] = ggl[k] - suml;
      if (std::abs(diL[j]) < 1e-10)
      {
        std::cerr << "Error: Near-zero diagonal element diL[" << j
                  << "] = " << diL[j] << std::endl;
        exit(1);
      }
      U[k] = (ggu[k] - sumu) / diL[j];
      sumDi += L[k] * U[k];
    }

    // Вычисляем diL[i] после обработки всех элементов строки
    diL[i] = di[i] - sumDi;
    if (std::abs(diL[i]) < 1e-10)
    {
      std::cerr << "Error: Near-zero diagonal element diL[" << i
                << "] = " << diL[i] << std::endl;
      exit(1);
    }
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
