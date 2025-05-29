// #include "funcs.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>
using namespace std;

/*структура узла сетки*/
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

/* ввод-вывод
--------------------------------------------------------------------*/
// ввод данных
void input(vector<nd> &mesh, vector<el> &elList, double &gamma, int &fnum)
{
  int n, node;
  double r, z, l;

  // количество узлов и сами узлы в файле nds.txt
  ifstream in("../data/nds.txt");
  in >> n;
  mesh.resize(n);
  for (int i = 0; i < n; i++)
  {
    in >> r >> z;
    mesh[i] = nd(i, r, z);
  }
  in.close();

  // количество элементов и сами элементы в файле els.txt
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

  // значения лямбда на каждом элементе в файле lambda.txt
  in.open("../data/lambda.txt");
  for (int i = 0; i < n; i++)
    in >> elList[i].lambda;
  in.close();

  // значение гамма в файле gamma.txt
  in.open("../data/gamma.txt");
  in >> gamma;
  in.close();

  // номер функции f в файле fnum.txt
  in.open("../data/fnum.txt");
  in >> fnum;
  in.close();
}

// вывод результатов
void Out(vector<double> &q, vector<nd> &mesh)
{
  cout << scientific << setprecision(6);
  for (int i = 0; i < q.size(); ++i)
  {
    cout << "node " << mesh[i].gl_num << " (r = " << mesh[i].r
         << ", z = " << mesh[i].z << "): q = " << q[i] << endl;
  }
}

/*функции для вычислений
-------------------------------------------------------------------*/

// вычисление модуля определителя матрицы D
double detD(el e)
{
  return abs((e.nds[1].r - e.nds[0].r) * (e.nds[2].z - e.nds[0].z) -
             (e.nds[2].r - e.nds[0].r) * (e.nds[1].z - e.nds[0].z));
}

// вычисление длины ребра (n1, n2)
double len(nd n1, nd n2)
{
  return (sqrt((n2.r - n1.r) * (n2.r - n1.r) + (n2.z - n1.z) * (n2.z - n1.z)));
}

// вычисление a!
long fact(int a)
{
  long f = 1;
  for (int i = 1; i <= a; i++)
    f *= i;
  return f;
}

// интегрирование произведения L1*L2*L3, где nv[] - массив степеней
// соответствующих функций
double intL(int nv[], double det)
{
  double num = 1, den = 2;
  for (int i = 0; i < 3; i++)
  {
    num *= fact(nv[i]);
    den += nv[i];
  }
  den = fact(den);
  return (num * det) / (den * 2);
}

// умножение матрицы на вектор
void mmxv(vector<vector<double>> &C, vector<double> &f, vector<double> &b)
{
  int n = 3;

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      b[i] += C[i][j] * f[j];
}

/*значения f и u в зависимости от параметра
--------------------------------------------------------*/

// возвращает функцию f в зависимости от параметра
double f(int s, nd n, el e)
{
  switch (s)
  {
  case 0:
    return 0;
  case 1:
    return 1;
  case 2:
    return n.r - 1 / n.r;
  case 3:
    return 10 * n.r + 4 * n.z; // n.z
    // return n.r - 1 / n.r + n.z;
  case 4:
    return n.r * n.r - 4;
  // return n.r * n.r * n.r - 6 * n.r;
  case 5:
    return n.z * n.z - 2;
    // case 6:
    // return n.z * n.r - n.z / n.r;
  case 7:
    return 2 * sin(n.z);
    // return (n.z * n.z * n.z) - 6 * n.z;
    // return (1 / (n.z * n.z)) + log(n.z);
    // return 2 * exp(n.r);
  }
}
// возвращает функцию u в зависимости от параметра
double u(int s, nd n)
{
  switch (s)
  {
  case 0:
    return 0;
  case 1:
    return 0.1; // 1
    // return n.r + n.z;
  case 2:
    return n.r;
  case 3:
    return 5 * n.r + 2;
  case 4:
    return n.r * n.r;
    // return n.r * n.r * n.r;
  case 5:
    return n.z * n.z;
    // case 6:
    // return n.r * n.z;
  case 7:
    return sin(n.z);
    // return exp(n.r);
    // return log(n.z);
    // return n.z * n.z * n.z;
  }
}

/* Содержит функции для вычисления локальных матриц жесткости и массы.
-------------------------------------------------------------------------*/

// вычисление элемента матрицы массы, возвращает M_ij
double Mij(int i, int j, el e, double det)
{
  int nv[3];
  double sum = 0;
  for (int m = 0; m < 3; m++)
  {
    nv[0] = 0;
    nv[1] = 0;
    nv[2] = 0;
    nv[i]++;
    nv[j]++;
    nv[m]++;
    sum += e.nds[m].r * intL(nv, det);
  }
  return sum;
}

// получение локальной матрицы масс M, записывает результат в передаваемую
// матрицу.
void getM(vector<vector<double>> &M, double gamma, el e)
{
  double det = detD(e);
  int i, j;

  for (i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      M[i][j] = gamma * Mij(i, j, e, det);
}

// получение локальной матрицы жесткости G, записывает результат в переданную
// матрицу
void getG(vector<vector<double>> &G, el e)
{
  vector<double> a1(3), a2(3);
  double det;
  a1[0] = (e.nds[1].z - e.nds[2].z); // z2 - z3
  a1[1] = (e.nds[2].z - e.nds[0].z); // z3 - z1
  a1[2] = (e.nds[0].z - e.nds[1].z); // z1 - z2
  a2[0] = (e.nds[2].r - e.nds[1].r); // r3 - r2
  a2[1] = (e.nds[0].r - e.nds[2].r); // r1 - r3
  a2[2] = (e.nds[1].r - e.nds[0].r); // r2 - r1

  det = (detD(e) * e.lambda) * (e.nds[0].r + e.nds[1].r + e.nds[2].r) / 12;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
      G[i][j] = det * (a1[i] * a1[j] + a2[i] * a2[j]);
  }
}

// получение вектора b, результат записывается в переданный вектор
void getb(vector<double> &b, el e, double gamma, int s)
{
  double sum = 0;
  double det = detD(e);

  for (int i = 0; i < 3; i++)
  {
    for (long k = 0; k < 3; k++)
      sum += f(s, e.nds[k], e) * Mij(i, k, e, det);
    b[i] = sum;
    sum = 0;
  }
}

/*основные функции для добавления в глобальную матрицу и решения
---------------------------------------------------------*/
// поэлементная сборка глобальной матрицы A
void addA(vector<vector<double>> &A, vector<vector<double>> l, el e)
{

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      A[e.nds[i].gl_num][e.nds[j].gl_num] += l[i][j];
}

void addB(vector<double> &B, vector<double> b, el e)
{
  for (int i = 0; i < 3; i++)
    B[e.nds[i].gl_num] += b[i];
}

// Функция итерационного метода Якоби
void slv(vector<vector<double>> &A, vector<double> &B, vector<double> &X,
         double tlrc, int maxit)
{
  int n = A.size();
  vector<double> X_prev(n);

  for (int it = 0; it < maxit; ++it)
  {
    X_prev = X;

    for (int i = 0; i < n; ++i)
    {
      double sum = 0.0;

      for (int j = 0; j < n; ++j)
        if (j != i)
          sum += A[i][j] * X_prev[j];

      X[i] = (B[i] - sum) / A[i][i];
    }

    double err = 0.0;
    for (int i = 0; i < n; ++i)
      err = max(err, fabs(X[i] - X_prev[i]));

    if (err < tlrc)
    {
      cout << it + 1 << "iters." << endl;
      break;
    }
  }
}

/*учет первых краевых условий
---------------------------------------------------------*/
// учет первых краевых условий. A - глобальная матрица, b - глобавльный вектор,
// mesh - вектор с узлами сетки
void c1(vector<vector<double>> &A, vector<double> &b, vector<nd> mesh)
{
  int p;
  ifstream in("../data/c1.txt");
  in >> p;
  for (int i = 0; i < p; i++)
  {
    int v1, v2, s;
    in >> v1 >> v2 >> s;
    A[v1][v1] = 1;
    A[v2][v2] = 1;
    b[v1] = u(s, mesh[v1]);
    b[v2] = u(s, mesh[v2]);

    for (int j = 0; j < A.size(); j++)
    {
      if (j != v1)
        A[v1][j] = 0;
      if (j != v2)
        A[v2][j] = 0;
    }
  }
  in.close();
}

/*основная функция
---------------------------------------------------------*/
int main()
{
  vector<nd> mesh;
  vector<el> elList;
  vector<vector<double>> M(3), G(3);
  double gamma = 0; // Константа
  int fnum = 0;     // Число граничных условий

  // Чтение сетки и начальных данных
  input(mesh, elList, gamma, fnum);
  int n = mesh.size(); // Размер сетки

  vector<double> locb(3, 0), B(n), q(n, 0);
  vector<vector<double>> A(n);

  // изменяем размер векторов
  for (int i = 0; i < 3; i++)
  {
    M[i].resize(3);
    G[i].resize(3);
  }
  for (int i = 0; i < A.size(); i++)
  {
    A[i].resize(n);
    for (int j = 0; j < A.size(); j++)
      A[i][j] = 0;
  }

  // основной процесс расчета локальных матриц и добаления в глобальные
  for (int i = 0; i < elList.size(); i++)
  {
    getM(M, gamma, elList[i]);
    addA(A, M, elList[i]);
    getG(G, elList[i]);
    addA(A, G, elList[i]);
    getb(locb, elList[i], gamma, fnum);
    addB(B, locb, elList[i]);
  }

  // учет первых краевых, решение, вывод решения
  c1(A, B, mesh);
  slv(A, B, q, 1e-16, 10000);
  Out(q, mesh);

  return 0;
}
