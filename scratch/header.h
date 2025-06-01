#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>
using namespace std;

/*============================= DATA STRUCTURES =============================*/
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

// time mesh structure
struct TimeMesh
{
  vector<double> t{};
  vector<double> tNew{};

  vector<vector<double>> q{};
  vector<double> qti{};
  vector<double> qti_1{};
  vector<double> qti_2{};

  void Swap()
  {
    std::swap(qti_2, qti_1);
    std::swap(qti_1, qti);
  };
  void SaveWeights(int ti) { q[ti] = qti; };
};

// time mesh funcs

// sparse matrix structure
struct SparseMatrix
{
  vector<int> ig{}, jg{};
  vector<double> ggl{}, ggu{}, di{};
};

// SLAE structure
struct SLAE
{
  SparseMatrix A{};
  vector<double> q{}, b{};
};

// LOS structure
struct LOS
{
  vector<double> r1{}, rk{}, z1{}, p1{}, Ar{}, p{}, mult{};
};

// Boundary conditions structure
struct BoundaryConditions
{
  int type{}, function{};
  vector<int> VerticesNumbers{};
};

/*============================== I/O  FUNCTIONS =============================*/

// Input data
void input(vector<nd> &mesh, vector<el> &elList, double &sigma, int &fnum);
// Read time mesh
void readTimeMesh(TimeMesh &time);
// Read time steps
void readSplitTimeMesh(TimeMesh &time);
// boundaries input
void readBoundaryCondition(vector<BoundaryConditions> &cond);

void PrintSolution(const TimeMesh &timemesh, const vector<nd> &mesh);

// Output results
// void output(const vector<double> &q, const vector<nd> &mesh, double t);

/*========================= LOCAL MATRICES FUNCTIONS =======================*/

// Compute determinant
double detD(el &e);
// Compute factorial
long fact(int a);
// Integrate L1*L2*L3
double intL(int nv[], double det);
// Compute mass matrix element M_ij
double Mij(int i, int j, el e, double det);
// Get local mass matrix M
void getM(vector<vector<double>> &M, double sigma, el &e);
// Get local stiffness matrix G
void getG(vector<vector<double>> &G, el e);
// Get local load vector b
void getb(vector<double> &b, el e, double t, double gamma, int flag);

/*======================== SPARSE MATRICES FUNCTIONS =======================*/
// Get portrait sparse matrix
void GetPortraitSparseMatrix(vector<nd> &mesh, vector<el> &elList, SLAE &slae);
// Get global matrix and vector
// проблема с добавлением лок в глобальную чет не так с этим вектором узлов
void addLocalMatrixToGlobal(SparseMatrix &A, el elem,
                            vector<vector<double>> &localMatrix);
void addLocalVectorToGlobal(vector<double> &b, el &elem,
                            vector<double> &bLocal);
void addElemToGlobalMatrix(SparseMatrix &A, int i, int j, double elem);
void multiplyMatrixToVector(vector<vector<double>> &matrix, vector<double> &vec,
                            vector<double> &result, el &elem);
void multiplyMatrixToCoef(vector<vector<double>> &matrix, double coef,
                          vector<vector<double>> &resultMatrix);
void multiplyVectorToCoef(vector<double> &vector, double coef);

/*============================== LOS FUNCTIONS ==============================*/

void calcLU(SLAE &slae, SLAE &LU);
void localOptimalSchemeLU(SLAE &slae, SLAE &LU, LOS &v, int maxIter,
                          double eps);
void calcY(SLAE &LU, vector<double> &b, vector<double> &y);
void calcX(SLAE &LU, vector<double> &y, vector<double> &x);
void multOfMatrix(SparseMatrix &A, vector<double> &x, vector<double> &F);
void calcDiscrepancy(SLAE &slae, LOS &v, vector<double> &x, double &normb);
void clearSLAE(SLAE &slae);

/*============================= SOLVER FUNCTIONS =============================*/

// solver
void Solve(vector<nd> &mesh, vector<el> &elList, TimeMesh &timemesh, SLAE &slae,
           vector<BoundaryConditions> &conds, int flag, double sigma);

void GetGlobalMatrixAndVector(vector<nd> &mesh, vector<el> elList,
                              TimeMesh &timemesh, SLAE &slae,
                              vector<BoundaryConditions> &cond, int ti,
                              int flag, double sigma);

// init q0 q1
void getWeightsInitU(TimeMesh &timemesh, vector<nd> &mesh);

void addFirstBoundaryCondition(SLAE &slae, vector<BoundaryConditions> &cond,
                               vector<nd> mesh, double tValue);

void stringMatrixInNull(SparseMatrix &A, int i);

/*============================== F/U FUNCTIONS =============================*/

double f(nd &n, double t, int flag);
double u(nd &n, double t, int flag);
double uInit(nd &n, double t, int flag);
void DebugSolver(vector<nd> &mesh, vector<el> &elList, TimeMesh &timemesh,
                 SLAE &slae, vector<BoundaryConditions> &conds, int flag,
                 double sigma);
// void PrintAllInfo(const SLAE &slae, const TimeMesh &tm, const LOS &los);

// // add switch!!!!!!!!1
// double uInit(nd &n, double t)
// { // примерная функция, уточнить
//   return n.r * n.r + n.z * n.z + t * t;
// }
// double sigma(nd &n, double t)
// {
//   // Примерная функция сигмы, нужно уточнить
//   return n.r * n.r + n.z * n.z + t * t;
// }
// double F(nd &n, double t)
// {
//   // Примерная функция, нужно уточнить
//   return n.r * n.r + n.z * n.z + t * t;
// }
