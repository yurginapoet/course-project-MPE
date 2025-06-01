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
    qti_2[i] = uInit(mesh[i], timemesh.tNew[0], flag);
    qti_1[i] = uInit(mesh[i], timemesh.tNew[1], flag);
  }
}

// solver
void Solve(vector<nd> &mesh, vector<el> &elList, TimeMesh &timemesh, SLAE &slae,
           vector<BoundaryConditions> &conds, int flag, double sigma)
{
  int tSize = timemesh.tNew.size();

  getWeightsInitU(timemesh, mesh, flag);
  timemesh.q[0] = timemesh.qti_2;
  timemesh.q[1] = timemesh.qti_1;

  for (int ti = 2; ti < tSize; ti++)
  {
    SLAE LU{};
    LOS vectors{};

    GetGlobalMatrixAndVector(mesh, elList, timemesh, slae, conds, ti, flag,
                             sigma);

    calcLU(slae, LU);
    localOptimalSchemeLU(slae, LU, vectors, 10000, 1e-14);

    timemesh.qti = slae.q;
    timemesh.SaveWeights(ti);
    timemesh.Swap();

    clearSLAE(slae);
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
    multiplyMatrixToCoef(M, (deltaT + deltaT0) / (deltaT * deltaT0),
                         tempMatrix);
    // ДОБАВЛЯЕМ ЛОКАЛЬНЫЕ МАТРИЦЫ В ГЛОБАЛЬНУЮ
    addLocalMatrixToGlobal(A, elem, tempMatrix);
    addLocalMatrixToGlobal(A, elem, G);

    // ДОБАВЛЯЕМ ЛОКАЛЬНЫЙ ВЕКТОР В ГЛОБАЛЬНЫЙ
    addLocalVectorToGlobal(b, elem, locb);
    fill(tempVector.begin(), tempVector.end(), 0.0);
    multiplyVectorToCoef(tempVector, -deltaT0 / (deltaT * deltaT1));
    multiplyMatrixToVector(M, timemesh.qti_2, tempVector, elem);

    addLocalVectorToGlobal(b, elem, tempVector);
    fill(tempVector.begin(), tempVector.end(), 0.0);
    multiplyMatrixToVector(M, timemesh.qti_1, tempVector, elem);
    multiplyVectorToCoef(tempVector, deltaT / (deltaT1 * deltaT0));
    addLocalVectorToGlobal(b, elem, tempVector);
  }

  // addSecondBoundaryCondition(slae, cond, vertexCoord, t);
  addFirstBoundaryCondition(slae, cond, mesh, t);
}

// void addFirstBoundaryCondition(SLAE &slae, vector<BoundaryConditions> &cond,
//                                vector<nd> mesh, double tValue)
// {
//   auto &A = slae.A;
//   auto &b = slae.b;

//   for (auto &condition : cond)
//     if (condition.type == 1)
//     {
//       auto &VerticesNumbers = condition.VerticesNumbers;
//       int flag = condition.function;

//       const int numVertex = VerticesNumbers.size();

//       for (int i = 0; i < numVertex; i++)
//       {
//         int globalNum = VerticesNumbers[i];
//         stringMatrixInNull(A, globalNum);
//         b[globalNum] = u(mesh[globalNum], tValue, flag);
//       }
//     }
// }

#include "header.h"

void addFirstBoundaryCondition(SLAE &slae, vector<BoundaryConditions> &cond,
                               vector<nd> mesh, double tValue)
{
  auto &A = slae.A;
  auto &b = slae.b;

  std::ofstream debugOut("debug_boundary.txt", std::ios::app);
  debugOut << std::fixed << std::setprecision(10);
  debugOut << "Applying Boundary Conditions at t = " << tValue << "\n";

  // Сохраняем копию вектора b до применения граничных условий
  std::vector<double> b_before = b;
  debugOut << "Vector b before boundary conditions:\n";
  for (size_t i = 0; i < b.size(); ++i)
    debugOut << "b[" << i << "] = " << b[i] << "\n";

  for (auto &condition : cond)
  {
    if (condition.type == 1)
    {
      auto &VerticesNumbers = condition.VerticesNumbers;
      int flag = condition.function;

      debugOut << "Boundary condition: type = " << condition.type
               << ", function = " << flag
               << ", numVertex = " << VerticesNumbers.size() << "\n";
      debugOut << "Vertices: ";
      for (int v : VerticesNumbers)
        debugOut << v << " ";
      debugOut << "\n";

      // Проверка, что не все узлы граничные
      if (VerticesNumbers.size() == mesh.size())
      {
        debugOut << "[WARNING] All nodes are boundary nodes! This may zero out "
                    "vector b.\n";
      }

      for (int i = 0; i < VerticesNumbers.size(); i++)
      {
        int globalNum = VerticesNumbers[i];
        if (globalNum < 0 || globalNum >= static_cast<int>(mesh.size()))
        {
          debugOut << "[ERROR] Invalid node index: " << globalNum << "\n";
          std::cerr << "Error: Invalid node index in boundary conditions: "
                    << globalNum << std::endl;
          exit(1);
        }

        double uValue = u(mesh[globalNum], tValue, flag);
        debugOut << "Node " << globalNum << " (r = " << mesh[globalNum].r
                 << ", z = " << mesh[globalNum].z << "): u = " << uValue;
        if (std::isnan(uValue) || std::isinf(uValue))
        {
          debugOut << " [ERROR: NaN or Inf in u]\n";
          std::cerr << "Error: u returns NaN or Inf for node " << globalNum
                    << std::endl;
          exit(1);
        }
        debugOut << "\n";

        // Обнуляем строку и столбец в матрице A и устанавливаем диагональный
        // элемент
        stringMatrixInNull(A, globalNum);
        // Устанавливаем значение в векторе b
        b[globalNum] = uValue;
      }
    }
  }

  debugOut << "Vector b after boundary conditions:\n";
  for (size_t i = 0; i < b.size(); ++i)
  {
    debugOut << "b[" << i << "] = " << b[i];
    if (b[i] == 0.0 && b_before[i] != 0.0)
      debugOut << " [WARNING: Value zeroed from " << b_before[i] << "]";
    debugOut << "\n";
  }

  // Проверка, что вектор b не стал полностью нулевым
  double normb = 0.0;
  for (double val : b)
    normb += val * val;
  if (normb < 1e-10)
  {
    debugOut << "[ERROR] Vector b is zero after boundary conditions (normb = "
             << normb << ")\n";
    std::cerr << "Error: Vector b is zero after applying boundary conditions\n";
    exit(1);
  }

  debugOut << "\n";
  debugOut.close();
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

#include "header.h"
#include <iomanip>

void DebugSolver(vector<nd> &mesh, vector<el> &elList, TimeMesh &timemesh,
                 SLAE &slae, vector<BoundaryConditions> &conds, int flag,
                 double sigma)
{
  // Открываем файл для записи отладочной информации
  ofstream debugOut("debug_solver_output.txt");
  debugOut << std::fixed << std::setprecision(10);

  debugOut << "=== Debugging Solver ===\n";
  debugOut << "Number of nodes: " << mesh.size() << "\n";
  debugOut << "Number of elements: " << elList.size() << "\n";
  debugOut << "Number of time steps: " << timemesh.tNew.size() << "\n";
  debugOut << "Sigma: " << sigma << "\n";
  debugOut << "Flag: " << flag << "\n\n";

  // Проверка начальных условий
  debugOut << "=== Initial Conditions ===\n";
  getWeightsInitU(timemesh, mesh, flag);
  debugOut << "qti_2 (t = " << timemesh.tNew[0] << "):\n";
  for (size_t i = 0; i < timemesh.qti_2.size(); ++i)
  {
    debugOut << "Node " << i << ": " << timemesh.qti_2[i];
    if (std::isnan(timemesh.qti_2[i]) || std::isinf(timemesh.qti_2[i]))
      debugOut << " [WARNING: NaN or Inf]";
    debugOut << "\n";
  }
  debugOut << "qti_1 (t = " << timemesh.tNew[1] << "):\n";
  for (size_t i = 0; i < timemesh.qti_1.size(); ++i)
  {
    debugOut << "Node " << i << ": " << timemesh.qti_1[i];
    if (std::isnan(timemesh.qti_1[i]) || std::isinf(timemesh.qti_1[i]))
      debugOut << " [WARNING: NaN or Inf]";
    debugOut << "\n";
  }
  timemesh.q[0] = timemesh.qti_2;
  timemesh.q[1] = timemesh.qti_1;

  // Основной цикл по времени
  int tSize = timemesh.tNew.size();
  for (int ti = 2; ti < tSize; ++ti)
  {
    debugOut << "\n=== Time Step " << ti << ", t = " << timemesh.tNew[ti]
             << " ===\n";

    // Проверка временных шагов
    double deltaT = timemesh.tNew[ti] - timemesh.tNew[ti - 2];
    double deltaT0 = timemesh.tNew[ti] - timemesh.tNew[ti - 1];
    double deltaT1 = timemesh.tNew[ti - 1] - timemesh.tNew[ti - 2];
    debugOut << "deltaT: " << deltaT << ", deltaT0: " << deltaT0
             << ", deltaT1: " << deltaT1 << "\n";
    if (std::abs(deltaT) < 1e-10 || std::abs(deltaT0) < 1e-10 ||
        std::abs(deltaT1) < 1e-10)
    {
      debugOut << "[ERROR] Time step is too small or zero!\n";
    }

    // Вызов GetGlobalMatrixAndVector
    debugOut << "\n=== Debugging GetGlobalMatrixAndVector ===\n";
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

    debugOut << "Processing elements:\n";
    for (size_t e = 0; e < elList.size(); ++e)
    {
      auto &elem = elList[e];
      debugOut << "Element " << e << " (nodes: " << elem.nds[0].gl_num << ", "
               << elem.nds[1].gl_num << ", " << elem.nds[2].gl_num << "):\n";

      // Вычисление локальных матриц и вектора
      getM(M, sigma, elem);
      getG(G, elem);
      getb(locb, elem, t, sigma, flag);

      // Проверка локальной матрицы масс M
      debugOut << "Local Mass Matrix M:\n";
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          debugOut << M[i][j] << " ";
          if (std::isnan(M[i][j]) || std::isinf(M[i][j]))
            debugOut << "[WARNING: NaN or Inf] ";
        }
        debugOut << "\n";
      }

      // Проверка локальной матрицы жесткости G
      debugOut << "Local Stiffness Matrix G:\n";
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          debugOut << G[i][j] << " ";
          if (std::isnan(G[i][j]) || std::isinf(G[i][j]))
            debugOut << "[WARNING: NaN or Inf] ";
        }
        debugOut << "\n";
      }

      // Проверка локального вектора b
      debugOut << "Local Load Vector b:\n";
      for (int i = 0; i < 3; ++i)
      {
        debugOut << locb[i] << " ";
        if (std::isnan(locb[i]) || std::isinf(locb[i]))
          debugOut << "[WARNING: NaN or Inf] ";
      }
      debugOut << "\n";

      // Умножение матрицы M на коэффициент
      double coef = (deltaT + deltaT0) / (deltaT * deltaT0);
      debugOut << "Coefficient for M: " << coef << "\n";
      if (std::isnan(coef) || std::isinf(coef))
        debugOut << "[ERROR] Invalid coefficient for M!\n";
      multiplyMatrixToCoef(M, coef, tempMatrix);

      // Добавление в глобальную матрицу
      addLocalMatrixToGlobal(A, elem, tempMatrix);
      addLocalMatrixToGlobal(A, elem, G);

      // Добавление локального вектора b
      addLocalVectorToGlobal(b, elem, locb);

      // Обработка qti_2 и qti_1
      fill(tempVector.begin(), tempVector.end(), 0.0);
      multiplyMatrixToVector(M, timemesh.qti_2, tempVector, elem);
      multiplyVectorToCoef(tempVector, -deltaT0 / (deltaT * deltaT1));
      debugOut << "tempVector after qti_2 contribution:\n";
      for (int i = 0; i < 3; ++i)
      {
        debugOut << tempVector[i] << " ";
        if (std::isnan(tempVector[i]) || std::isinf(tempVector[i]))
          debugOut << "[WARNING: NaN or Inf] ";
      }
      debugOut << "\n";
      addLocalVectorToGlobal(b, elem, tempVector);

      fill(tempVector.begin(), tempVector.end(), 0.0);
      multiplyMatrixToVector(M, timemesh.qti_1, tempVector, elem);
      multiplyVectorToCoef(tempVector, deltaT / (deltaT1 * deltaT0));
      debugOut << "tempVector after qti_1 contribution:\n";
      for (int i = 0; i < 3; ++i)
      {
        debugOut << tempVector[i] << " ";
        if (std::isnan(tempVector[i]) || std::isinf(tempVector[i]))
          debugOut << "[WARNING: NaN or Inf] ";
      }
      debugOut << "\n";
      addLocalVectorToGlobal(b, elem, tempVector);
    }

    // Вывод глобальной матрицы и вектора
    debugOut << "\nGlobal Matrix A (diagonal and non-zero elements):\n";
    for (size_t i = 0; i < A.di.size(); ++i)
    {
      debugOut << "Row " << i << ": di = " << A.di[i];
      if (std::isnan(A.di[i]) || std::isinf(A.di[i]))
        debugOut << " [WARNING: NaN or Inf]";
      debugOut << "\n";
      int i0 = A.ig[i], i1 = A.ig[i + 1];
      for (int k = i0; k < i1; ++k)
      {
        debugOut << "  jg[" << k << "] = " << A.jg[k] << ", ggl = " << A.ggl[k]
                 << ", ggu = " << A.ggu[k];
        if (std::isnan(A.ggl[k]) || std::isinf(A.ggl[k]) ||
            std::isnan(A.ggu[k]) || std::isinf(A.ggu[k]))
          debugOut << " [WARNING: NaN or Inf]";
        debugOut << "\n";
      }
    }

    debugOut << "Global Vector b:\n";
    for (size_t i = 0; i < b.size(); ++i)
    {
      debugOut << "b[" << i << "] = " << b[i];
      if (std::isnan(b[i]) || std::isinf(b[i]))
        debugOut << " [WARNING: NaN or Inf]";
      debugOut << "\n";
    }

    // Применение граничных условий
    debugOut << "\nApplying Boundary Conditions:\n";
    addFirstBoundaryCondition(slae, conds, mesh, t);

    // Проверка LU-разложения
    debugOut << "\n=== LU Decomposition ===\n";
    calcLU(slae, LU);
    for (size_t i = 0; i < LU.A.di.size(); ++i)
    {
      debugOut << "diL[" << i << "] = " << LU.A.di[i];
      if (std::abs(LU.A.di[i]) < 1e-10)
        debugOut << " [WARNING: Near-zero diagonal element]";
      if (std::isnan(LU.A.di[i]) || std::isinf(LU.A.di[i]))
        debugOut << " [WARNING: NaN or Inf]";
      debugOut << "\n";
    }

    // Решение системы методом LOS
    debugOut << "\n=== Solving with LOS ===\n";
    localOptimalSchemeLU(slae, LU, vectors, 10000, 1e-14);

    // Проверка решения
    debugOut << "Solution q:\n";
    for (size_t i = 0; i < slae.q.size(); ++i)
    {
      debugOut << "q[" << i << "] = " << slae.q[i];
      if (std::isnan(slae.q[i]) || std::isinf(slae.q[i]))
        debugOut << " [WARNING: NaN or Inf]";
      debugOut << "\n";
    }

    timemesh.qti = slae.q;
    timemesh.SaveWeights(ti);
    timemesh.Swap();
    clearSLAE(slae);
  }

  debugOut.close();
}