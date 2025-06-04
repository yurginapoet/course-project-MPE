#include "header.h"

/*=======================================================================*/

// вычисление модуля определителя матрицы D
double detD(el &e)
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

/*=======================================================================*/

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
void getM(vector<vector<double>> &M, double sigma, el &e)
{
  double det = detD(e);
  int i, j;

  for (i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      M[i][j] = sigma * Mij(i, j, e, det);
  // ???? точно ли тут должна быть сигма?
}

// получение локальной матрицы жесткости G, записывает результат в переданную
// матрицу
// что-то не так
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
// нужно переписать построение локального вектора, потому что там все
// по-другому??.
void getb(vector<double> &b, el e, double t, double gamma, int flag)
{
  double sum = 0;
  double det = detD(e);

  for (int i = 0; i < 3; i++)
  {
    for (long k = 0; k < 3; k++)
      sum += f(e.nds[k], t, flag) * Mij(i, k, e, det);
    b[i] = sum;
    sum = 0;
  }
}
