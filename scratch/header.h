#pragma once

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

// структура разреженной матрицы (не используется)
struct SparseMatrix
{
  vector<int> ia;    // Индексы строк (начало строк)
  vector<int> ja;    // Индексы столбцов
  vector<double> au; // Верхняя часть матрицы (значения на верхней диагонали)
  vector<double> al; // Нижняя часть матрицы (значения на нижней диагонали)
  vector<double> di; // Диагональные элементы (если они есть)
  int n;             // Размерность матрицы

  // Конструктор
  SparseMatrix(int size)
  {
    n = size;
    ia.resize(n + 1, 0);
    ja.resize(n, 0);
    au.resize(n, 0.0);
    al.resize(n, 0.0);
    di.resize(n, 0.0);
  }
};

/*========================= I/O  FUNCTIONS =========================*/

// Input data
void input(vector<nd> &mesh, vector<el> &elList, double &gamma, int &fnum,
           vector<double> &time);
// Output results
void output(const vector<double> &q, const vector<nd> &mesh, double t);

/*==================== LOCAL MATRICES FUNCTIONS ====================*/

double detD(const el &e);          // Compute determinant
long fact(int a);                  // Factorial
double intL(int nv[], double det); // Integrate L1*L2*L3
double Mij(int i, int j, const el &e,
           double det); // Compute mass matrix element M_ij
void getM(vector<vector<double>> &M, double gamma,
          const el &e);                            // Get local mass matrix
void getG(vector<vector<double>> &G, const el &e); // Get local stiffness matrix
void getb(vector<double> &b, const el &e, double gamma, int s, double t,
          const vector<double> &q1, const vector<double> &q2,
          double dt); // Get local load vector

/*=================== SPARSE MATRICES FUNCTIONS ===================*/

void buildPortrait(vector<vector<int>> &list, vector<int> &ig, vector<int> &jg,
                   const vector<el> &elList); // Build sparse matrix portrait

void addElementToGlobal(int i, int j, double elem, vector<double> &di,
                        vector<double> &gg, const vector<int> &ig,
                        const vector<int> &jg); // Add element to sparse matrix

void mult(const vector<double> &di, const vector<double> &gg,
          const vector<int> &ig, const vector<int> &jg, const vector<double> &x,
          vector<double> &y); // Matrix-vector multiplication for sparse matrix

void conjugateGradient(const vector<double> &di, const vector<double> &gg,
                       const vector<int> &ig, const vector<int> &jg,
                       const vector<double> &F, vector<double> &q, double eps,
                       int maxIter); // Conjugate gradient solver

/*========================= F/U FUNCTIONS =========================*/

// Function f (modify as needed)
double f(int s, const nd &n, const el &e)
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
    return 10 * n.r + 4 * n.z;
  case 4:
    return n.r * n.r - 4;
  case 5:
    return n.z * n.z - 2;
  case 7:
    return 2 * sin(n.z);
  default:
    return 0;
  }
}

// Function u for boundary conditions
double u(int s, const nd &n)
{
  switch (s)
  {
  case 0:
    return 0;
  case 1:
    return 0.1;
  case 2:
    return n.r;
  case 3:
    return 5 * n.r + 2;
  case 4:
    return n.r * n.r;
  case 5:
    return n.z * n.z;
  case 7:
    return sin(n.z);
  default:
    return 0;
  }
}
