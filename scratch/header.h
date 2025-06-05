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
// node structure
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

// element structure
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
  void save_ti(int ti) { q[ti] = qti; };
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

// boundary conditions structure
struct BoundaryConditions
{
  int type{}, function{};
  vector<int> VerticesNumbers{};
};

/*============================== I/O  FUNCTIONS =============================*/

// Input data
void input(vector<nd> &mesh, vector<el> &elList, double &sigma, int &fnum);
// Read time mesh
void input_timemesh(TimeMesh &time);
// Read time steps
void input_split_timemesh(TimeMesh &time);
// boundaries input
void input_boundary(vector<BoundaryConditions> &cond);

void print_solution(const TimeMesh &timemesh, const vector<nd> &mesh);

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
void portrait(vector<nd> &mesh, vector<el> &elList, SLAE &slae);
// Get global matrix and vector
void add_mx(SparseMatrix &A, el elem, vector<vector<double>> &localMatrix);
void add_vec(vector<double> &b, el &elem, vector<double> &bLocal);
void add_el(SparseMatrix &A, int i, int j, double elem);
void mult_mx_vec(vector<vector<double>> &matrix, vector<double> &vec,
                 vector<double> &result, el &elem);
void mult_mx_num(vector<vector<double>> &matrix, double coef,
                 vector<vector<double>> &resultMatrix);
void mult_vec_num(vector<double> &vector, double coef);

/*============================== LOS FUNCTIONS ==============================*/

void calcLU(SLAE &slae, SLAE &LU);
void losLU(SLAE &slae, SLAE &LU, LOS &v, int maxIter, double eps);
void calcY(SLAE &LU, vector<double> &b, vector<double> &y);
void calcX(SLAE &LU, vector<double> &y, vector<double> &x);
void mult_smx_vec(SparseMatrix &A, vector<double> &x, vector<double> &F);
void calc_discrepancy(SLAE &slae, LOS &v, vector<double> &x, double &normb);
void clearSLAE(SLAE &slae);

/*============================= SOLVER FUNCTIONS =============================*/

// solver
void Solver(vector<nd> &mesh, vector<el> &elList, TimeMesh &timemesh,
            SLAE &slae, vector<BoundaryConditions> &conds, int flag,
            double sigma);

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