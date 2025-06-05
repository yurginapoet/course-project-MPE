#include "header.h"

void getWeightsInitU(TimeMesh &timemesh, vector<nd> &mesh, int flag)
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
    qti_2[i] = u(mesh[i], timemesh.tNew[0], flag);
    qti_1[i] = u(mesh[i], timemesh.tNew[1], flag);
  }
}

void GetGlobalMatrixAndVector(vector<nd> &mesh, vector<el> elList,
                              TimeMesh &timemesh, SLAE &slae,
                              vector<BoundaryConditions> &cond, int ti,
                              int flag, double sigma)
{
  auto &A = slae.A;
  auto &b = slae.b;

  double t = timemesh.tNew[ti];
  double deltaT = timemesh.tNew[ti] - timemesh.tNew[ti - 2];
  double deltaT0 = timemesh.tNew[ti] - timemesh.tNew[ti - 1];
  double deltaT1 = timemesh.tNew[ti - 1] - timemesh.tNew[ti - 2];
  int sizeSlae = mesh.size();

  vector<vector<double>> M(3), G(3), tempMatrix(3);
  vector<double> locb(3, 0), tempVector(3, 0);

  // изменяем размер векторов
  for (int i = 0; i < 3; i++)
  {
    M[i].resize(3);
    G[i].resize(3);
    tempMatrix[i].resize(3);
  }

  for (auto &elem : elList)
  {
    // НУЖНО ДОБАВИТЬ СИГМУ

    // ПРОСТО ТУПО ВЫЧИСЛЯЕМ ЛОКАЛЬНЫЕ МАТРИЦЫ И ВЕКТОРЫ
    getG(G, elem);
    getM(M, sigma, elem);
    getb(locb, elem, t, sigma, flag);

    // Умножаем локальные матрицы на коэффициенты
    mult_mx_num(M, (deltaT + deltaT0) / (deltaT * deltaT0), tempMatrix);
    // ДОБАВЛЯЕМ ЛОКАЛЬНЫЕ МАТРИЦЫ В ГЛОБАЛЬНУЮ
    add_mx(A, elem, tempMatrix);
    add_mx(A, elem, G);

    // ДОБАВЛЯЕМ ЛОКАЛЬНЫЙ ВЕКТОР В ГЛОБАЛЬНЫЙ
    add_vec(b, elem, locb);
    fill(tempVector.begin(), tempVector.end(), 0.0);
    mult_vec_num(tempVector, -deltaT0 / (deltaT * deltaT1));
    mult_mx_vec(M, timemesh.qti_2, tempVector, elem);

    add_vec(b, elem, tempVector);
    fill(tempVector.begin(), tempVector.end(), 0.0);
    mult_mx_vec(M, timemesh.qti_1, tempVector, elem);
    mult_vec_num(tempVector, deltaT / (deltaT1 * deltaT0));
    add_vec(b, elem, tempVector);
  }

  // addSecondBoundaryCondition(slae, cond, vertexCoord, t);
  addFirstBoundaryCondition(slae, cond, mesh, t);
}

void addFirstBoundaryCondition(SLAE &slae, vector<BoundaryConditions> &cond,
                               vector<nd> mesh, double tValue)
{
  auto &A = slae.A;
  auto &b = slae.b;

  for (auto &condition : cond)
    if (condition.type == 1)
    {
      auto &VerticesNumbers = condition.VerticesNumbers;
      int flag = condition.function;

      const int numVertex = VerticesNumbers.size();

      for (int i = 0; i < numVertex; i++)
      {
        int globalNum = VerticesNumbers[i];
        stringMatrixInNull(A, globalNum);
        b[globalNum] = u(mesh[globalNum], tValue, flag);
      }
    }
}

void stringMatrixInNull(SparseMatrix &A, int i)
{
  auto &ig = A.ig, &jg = A.jg;
  auto &ggl = A.ggl, &ggu = A.ggu, &di = A.di;
  const int size = ig[di.size()];

  int i0 = ig[i], i1 = ig[i + 1];

  for (i0; i0 < i1; i0++)
    ggl[i0] = 0.;

  int j0 = i1;

  for (j0; j0 < size; j0++)
    if (jg[j0] == i)
      ggu[j0] = 0.;

  di[i] = 1.;
}

void Solver(vector<nd> &mesh, vector<el> &elList, TimeMesh &timemesh,
            SLAE &slae, vector<BoundaryConditions> &conds, int flag,
            double sigma)
{
  getWeightsInitU(timemesh, mesh, flag);

  timemesh.q[0] = timemesh.qti_2;
  timemesh.q[1] = timemesh.qti_1;

  int tSize = timemesh.tNew.size();
  for (int ti = 2; ti < tSize; ++ti)
  {
    double deltaT = timemesh.tNew[ti] - timemesh.tNew[ti - 2];
    double deltaT0 = timemesh.tNew[ti] - timemesh.tNew[ti - 1];
    double deltaT1 = timemesh.tNew[ti - 1] - timemesh.tNew[ti - 2];

    SLAE LU{};
    LOS vectors{};
    auto &A = slae.A;
    auto &b = slae.b;
    double t = timemesh.tNew[ti];
    int sizeSlae = mesh.size();

    vector<vector<double>> M(3, vector<double>(3, 0.0));
    vector<vector<double>> G(3, vector<double>(3, 0.0));
    vector<vector<double>> tempMatrix(3, vector<double>(3, 0.0));
    vector<double> locb(3, 0.0), tempVector(3, 0.0);

    for (size_t e = 0; e < elList.size(); ++e)
    {
      auto &elem = elList[e];
      getM(M, sigma, elem);
      getG(G, elem);
      getb(locb, elem, t, sigma, flag);

      double coef = (deltaT + deltaT0) / (deltaT * deltaT0);
      mult_mx_num(M, coef, tempMatrix);
      add_mx(A, elem, tempMatrix);
      add_mx(A, elem, G);
      add_vec(b, elem, locb);

      fill(tempVector.begin(), tempVector.end(), 0.0);
      mult_mx_vec(M, timemesh.qti_2, tempVector, elem);
      mult_vec_num(tempVector, -deltaT0 / (deltaT * deltaT1));
      add_vec(b, elem, tempVector);

      fill(tempVector.begin(), tempVector.end(), 0.0);
      mult_mx_vec(M, timemesh.qti_1, tempVector, elem);
      mult_vec_num(tempVector, deltaT / (deltaT1 * deltaT0));
      add_vec(b, elem, tempVector);
    }

    addFirstBoundaryCondition(slae, conds, mesh, t);
    calcLU(slae, LU);
    losLU(slae, LU, vectors, 10000, 1e-30);

    timemesh.qti = slae.q;
    timemesh.save_ti(ti);
    timemesh.Swap();
    clearSLAE(slae);
  }
}