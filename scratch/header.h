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

// sparse matrix structure
struct smx
{
  vector<int> ig{}, jg{};
  vector<double> ggl{}, ggu{}, di{};
};

// SLAE structure
struct SLAE
{
  smx A{};
  vector<double> q{}, b{};
};

// LOS structure
struct LOS
{
  vector<double> r1{}, rk{}, z1{}, p1{}, Ar{}, p{}, mult{};
};

// boundary conditions structure
struct bc
{
  int type{}, function{};
  vector<int> ndnum{}; // number of nodes on the boundary
};

/*============================== I/O  FUNCTIONS =============================*/

// Input data
void input(vector<nd> &mesh, vector<el> &elList, double &sigma, int &fnum);
// Read time mesh
void input_timemesh(TimeMesh &time);
// Read time steps
void input_split_timemesh(TimeMesh &time);
// boundaries input
void input_boundary(vector<bc> &cond);
// print solution to the terminal
void print_solution(const TimeMesh &timemesh, const vector<nd> &mesh);
// output solution to the file
void out_solution(TimeMesh &timemesh, vector<nd> &mesh, int flag);

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
// Get portrait of the sparse matrix
void portrait(vector<nd> &mesh, vector<el> &elList, SLAE &slae);
// add matrix to the global sparse matrix
void add_mx(smx &A, el elem, vector<vector<double>> &localMatrix);
// add vector to the global vector
void add_vec(vector<double> &b, el &elem, vector<double> &bLocal);
// add element to the sparse matrix
void add_el(smx &A, int i, int j, double elem);
// multiplication of matrix and vector
void mult_mx_vec(vector<vector<double>> &matrix, vector<double> &vec,
                 vector<double> &result, el &elem);
// multiplication of matrix and number
void mult_mx_num(vector<vector<double>> &matrix, double coef,
                 vector<vector<double>> &resultMatrix);
// multiplication of vector and numver
void mult_vec_num(vector<double> &vector, double coef);

/*============================== LOS FUNCTIONS ==============================*/
// LU factorization
void calcLU(SLAE &slae, SLAE &LU);
// solve SLAE using LOS with LU decomposition
void losLU(SLAE &slae, SLAE &LU, LOS &v, int maxIter, double eps);
// forward substitution
void calcY(SLAE &LU, vector<double> &b, vector<double> &y);
// backward substitution
void calcX(SLAE &LU, vector<double> &y, vector<double> &x);

// multiplication of sparse matrix and vector
void mult_smx_vec(smx &A, vector<double> &x, vector<double> &F);
// calculate discrepancy
void calc_discrepancy(SLAE &slae, LOS &v, vector<double> &x, double &normb);
void clearSLAE(SLAE &slae);

/*============================= SOLVER FUNCTIONS =============================*/

// solver
void Solver(vector<nd> &mesh, vector<el> &elList, TimeMesh &timemesh,
            SLAE &slae, vector<bc> &conds, int flag, double sigma);

// get global sparse matrix
void get_global(vector<nd> &mesh, vector<el> elList, TimeMesh &timemesh,
                SLAE &slae, vector<bc> &cond, int ti, int flag, double sigma);

// initialize qti_2 and qti_1 for the first two time steps
void initq0q1(TimeMesh &timemesh, vector<nd> &mesh, int flag);
// apply first boundary conditions
void bc1(SLAE &slae, vector<bc> &cond, vector<nd> mesh, double tValue);
// clear line in sparse matrix leaving only diagonal element
void mx_clearline(smx &A, int i);

/*============================== F/U FUNCTIONS =============================*/

double f(nd &n, double t, int flag);
double u(nd &n, double t, int flag);