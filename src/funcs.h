#pragma once
#ifndef funcs_h
#define funcs_h

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

void input(vector<nd> &mesh, vector<el> &elList, double &gamma, int &fnum);
void Out(vector<double> &q, vector<nd> &mesh);

double detD(el e);
double len(nd n1, nd n2);
long fact(int a);
double intL(int nv[], double det);
void mmxv(vector<vector<double>> &C, vector<double> &f, vector<double> &b);

void getM(vector<vector<double>> &M, double gamma, el e);
void getG(vector<vector<double>> &G, el e);
void getb(vector<double> &b, el e, double gamma, int s);

double f(int s, nd n, el e);
double u(int s, nd n);

void addA(vector<vector<double>> &A, vector<vector<double>> l, el e);
void addB(vector<double> &B, vector<double> b, el e);
void slv(vector<vector<double>> &A, vector<double> &B, vector<double> &X,
         double tlrc, int maxit);
void c1(vector<vector<double>> &A, vector<double> &b, vector<nd> mesh);

#endif funcs_h