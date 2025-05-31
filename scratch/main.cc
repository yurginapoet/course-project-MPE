// #include "header.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>
using namespace std;

typedef std::vector<double> mvec;
// Vector addition
inline mvec operator+(mvec &a, const mvec &b)
{
  mvec res = a;
  for (int i = 0; i < a.size(); i++)
    res[i] += b[i];
  return res;
}

// Vector subtraction
inline mvec operator-(mvec &a, const mvec &b)
{
  mvec res = a;
  for (int i = 0; i < a.size(); i++)
    res[i] -= b[i];
  return res;
}

// Dot product
inline double operator*(const mvec &a, const mvec &b)
{
  double scalar = 0.0;
  for (int i = 0; i < a.size(); i++)
    scalar += a[i] * b[i];
  return scalar;
}

// Scalar-vector multiplication
inline mvec operator*(double c, mvec &a)
{
  mvec res = a;
  for (int i = 0; i < a.size(); i++)
    res[i] *= c;
  return res;
}

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

struct timelayer
{
  double t;      // value of current timestep
  int intervals; // number of intervals on current timestep
  double coef;   // coefficient of intervals on current time step
};

// sparse matrix
struct sparsematrix
{
  vector<int> ia;    // Индексы строк (начало строк)
  vector<int> ja;    // Индексы столбцов
  vector<double> au; // Верхняя часть матрицы (значения на верхней диагонали)
  vector<double> al; // Нижняя часть матрицы (значения на нижней диагонали)
  vector<double> di; // Диагональные элементы (если они есть)
  int n;             // Размерность матрицы

  // Конструктор
  sparsematrix(int size)
  {
    n = size;
    ia.resize(n + 1, 0);
    ja.resize(n, 0);
    au.resize(n, 0.0);
    al.resize(n, 0.0);
    di.resize(n, 0.0);
  }
};

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

void buildPortrait(vector<el> &elList, vector<nd> &mesh, sparsematrix &matrix)
{
  vector<vector<int>> list(mesh.size());
  list[0].push_back(0);

  // Go through all finite elements
  for (int ielem = 0; ielem < elList.size(); ielem++)
  {
    // Process pairs of nodes in the element
    for (int i = 0; i < 3; i++) // Each element has 3 nodes
    {
      for (int j = i + 1; j < 3; j++)
      {
        // Get global node numbers
        int node_i = elList[ielem].nds[i].gl_num;
        int node_j = elList[ielem].nds[j].gl_num;

        // Ensure node_i is smaller than node_j for upper triangle
        int insertPos = max(node_i, node_j);
        int element = min(node_i, node_j);

        bool isIn = false;

        // Check if element is already in the list
        for (int k = 0; k < list[insertPos].size() && !isIn; k++)
          if (element == list[insertPos][k])
            isIn = true;

        // If not found, add to the list
        if (!isIn)
          list[insertPos].push_back(element);
      }
    }
  }

  // Sort all lists in ascending order
  for (int i = 0; i < mesh.size(); i++)
    sort(list[i].begin(), list[i].end());

  // Form the ia array (row pointers)
  matrix.ia[0] = 0;
  matrix.ia[1] = 0;
  for (int i = 1; i < mesh.size(); i++)
    matrix.ia[i + 1] = matrix.ia[i] + list[i].size();

  // Form the ja array (column indices)
  matrix.ja.resize(matrix.ia[mesh.size()]);
  for (int i = 1, j = 0; i < mesh.size(); i++)
  {
    for (int k = 0; k < list[i].size(); k++)
      matrix.ja[j++] = list[i][k];
  }

  // Resize au, al, and di arrays to match ja
  matrix.au.resize(matrix.ja.size(), 0.0);
  matrix.al.resize(matrix.ja.size(), 0.0);
  matrix.di.resize(matrix.n, 0.0);
}

void toDense(const sparsematrix &matrix, const string &filename)
{
  vector<vector<double>> dense_matrix(matrix.n, vector<double>(matrix.n, 0.0));

  for (int i = 0; i < matrix.n; i++)
  {
    dense_matrix[i][i] = matrix.di[i];
    for (int j = matrix.ia[i]; j < matrix.ia[i + 1]; j++)
    {
      dense_matrix[i][matrix.ja[j]] = matrix.au[j];
      dense_matrix[matrix.ja[j]][i] = matrix.al[j];
    }
  }

  ofstream dense(filename);
  dense.precision(5);
  if (dense.is_open())
  {
    for (int i = 0; i < matrix.n; i++)
    {
      for (int j = 0; j <= i; j++)
        dense << left << setw(10) << dense_matrix[i][j];
      dense << endl << endl;
    }
  }
}

// Matrix-vector multiplication
void mult(const sparsematrix &matrix, const mvec &x, mvec &y)
{
  for (int i = 0; i < y.size(); i++)
    y[i] = 0;

  for (int i = 0; i < matrix.n; i++)
  {
    y[i] = matrix.di[i] * x[i];
    for (int k = matrix.ia[i]; k < matrix.ia[i + 1]; k++)
    {
      int j = matrix.ja[k];
      y[i] += matrix.au[k] * x[j];
      y[j] += matrix.al[k] * x[i];
    }
  }
}

// Euclidean norm
double EuclideanNorm(const mvec &x)
{
  double scalar = 0;
  for (int i = 0; i < x.size(); i++)
    scalar += x[i] * x[i];
  return sqrt(scalar);
}

// Conjugate gradient method
void CGM(const sparsematrix &matrix, mvec &q, mvec &F, double eps, int maxIter)
{
  int n = matrix.n;
  mvec um(n, 0.0), z(n, 0.0), r(n, 0.0);

  mult(matrix, q, um);
  r = F - um;
  z = r;

  double bnorm = EuclideanNorm(F);
  double residual = EuclideanNorm(r) / bnorm;

  if (residual > eps)
  {
    double scal1 = r * r;
    mult(matrix, z, um);
    double scal2 = um * z;
    double alfa = scal1 / scal2;

    q = q + alfa * z;
    r = r - alfa * um;

    double scal3 = r * r;
    double beta = scal3 / scal1;
    z = r + beta * z;
    residual = EuclideanNorm(r) / bnorm;
  }

  for (int k = 1; k < maxIter && residual > eps; k++)
  {
    double scal1 = r * r;
    mult(matrix, z, um);
    double scal2 = um * z;
    double alfa = scal1 / scal2;

    q = q + alfa * z;
    r = r - alfa * um;

    double scal3 = r * r;
    double beta = scal3 / scal1;
    z = r + beta * z;
    residual = EuclideanNorm(r) / bnorm;
  }
}

int main()
{
  // initialization
  vector<nd> mesh;
  vector<timelayer> time;
  vector<el> elList;
  double gamma;
  int fnum;

  // input
  input(mesh, elList, gamma, fnum, time);
  sparsematrix A(mesh.size());
  // building portrait
  buildPortrait(elList, mesh, A);

  return 0;
}