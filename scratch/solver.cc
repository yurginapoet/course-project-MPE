#include "header.h"

// solver
void Solve(vector<nd> &mesh, vector<el> &elList, TimeMesh &timemesh, SLAE &slae,
           vector<BoundaryConditions> &conds)
{
  int tSize = timemesh.tNew.size();

  getWeightsInitU(timemesh, mesh);
  timemesh.q[0] = timemesh.qti_2;
  timemesh.q[1] = timemesh.qti_1;

  for (int ti = 2; ti < tSize; ti++)
  {
    SLAE LU{};
    LOS vectors{};

    GetGlobalMatrixAndVector(mesh, elList, timemesh, slae, conds, ti);

    calcLU(slae, LU);
    localOptimalSchemeLU(slae, LU, vectors, 10000, 1e-14);

    timemesh.qti = slae.q;
    timemesh.saveWeights(ti);
    timemesh.swap();

    clearSLAE(slae);
  }
}

void GetGlobalMatrixAndVector(vector<nd> &mesh, vector<el> elList,
                              TimeMesh &timemesh, SLAE &slae,
                              vector<BoundaryConditions> &cond, int ti)
{
  auto &A = slae.A;
  auto &b = slae.b;

  double t = timemesh.tNew[ti];
  double deltaT = timemesh.tNew[ti] - timemesh.tNew[ti - 2];
  double deltaT0 = timemesh.tNew[ti] - timemesh.tNew[ti - 1];
  double deltaT1 = timemesh.tNew[ti - 1] - timemesh.tNew[ti - 2];
  int sizeSlae = mesh.size();

  vector<vector<double>> M(3), G(3), tempMatrix(3);
  vector<double> locb(3, 0), tempVector(3, 0), sigmas(3, 0);

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
    getSigmaWeights(mesh, elem, t);

    // ПРОСТО ТУПО ВЫЧИСЛЯЕМ ЛОКАЛЬНЫЕ МАТРИЦЫ И ВЕКТОРЫ
    getG(G, elem);
    getM(M, sigmas, elem);
    getb(locb, elem, t);

    // Умножаем локальные матрицы на коэффициенты
    multiplyMatrixToCoef(M, (deltaT + deltaT0) / (deltaT * deltaT0),
                         tempMatrix);
    // ДОБАВЛЯЕМ ЛОКАЛЬНЫЕ МАТРИЦЫ В ГЛОБАЛЬНУЮ
    addLocalMatrixToGlobal(A, elem.localVertex, tempMatrix);
    addLocalMatrixToGlobal(A, elem.localVertex, G);
    // ДОБАВЛЯЕМ ЛОКАЛЬНЫЙ ВЕКТОР В ГЛОБАЛЬНЫЙ
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