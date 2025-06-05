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

// Represents a node in the mesh.
// Stores the global node number and its coordinates (r, z) in a 2D axisymmetric
// domain.
struct nd
{
  int gl_num;  // Global node number
  double r, z; // Radial (r) and axial (z) coordinates
  nd(int _num, double _r, double _z)
      : gl_num(_num), r(_r), z(_z) {} // Constructor with initialization
  nd() = default;                     // Default constructor
};

// Represents a triangular element in the finite element mesh.
// Contains a material property (lambda) and an array of three nodes.
struct el
{
  double lambda; // Material property (e.g., thermal conductivity)
  nd nds[3];     // Array of three nodes forming the triangular element
};

// Represents the time mesh and solution data for time-dependent problems.
// Manages time steps and solution vectors across time layers.
struct TMesh
{
  vector<double> t{};         // Coarse time steps
  vector<double> tNew{};      // Refined time steps (after splitting)
  vector<vector<double>> q{}; // Solution at all time steps
  vector<double> qti{};       // Solution at current time step
  vector<double> qti_1{};     // Solution at previous time step
  vector<double> qti_2{};     // Solution at two time steps back

  // Swaps solution vectors to update time steps (qti_2 <- qti_1, qti_1 <- qti).
  void Swap()
  {
    std::swap(qti_2, qti_1);
    std::swap(qti_1, qti);
  };
  // Saves the current solution (qti) to the solution array q at time index ti.
  void save_ti(int ti) { q[ti] = qti; };
};

// Represents a sparse matrix in compressed storage format.
// Stores diagonal, lower triangular, and upper triangular elements separately.
struct smx
{
  vector<int> ig{};     // index of first non-zero element in each row
  vector<int> jg{};     // Column indices of non-zero elements
  vector<double> ggl{}; // Lower triangular elements
  vector<double> ggu{}; // Upper triangular elements
  vector<double> di{};  // Diagonal elements
};

// Represents a system of linear algebraic equations (SLAE).
// Contains the sparse matrix A, solution vector q, and right side vector b.
struct SLAE
{
  smx A{};            // Sparse matrix
  vector<double> q{}; // Solution vector
  vector<double> b{}; // Right- side vector
};

// Stores vectors used in the LOS solver.
struct LOS
{
  vector<double> r1{};   // Residual vector
  vector<double> rk{};   // Intermediate vector for iterations
  vector<double> z1{};   // Preconditioned residual
  vector<double> p1{};   // Search direction vector
  vector<double> Ar{};   // Matrix-vector product
  vector<double> p{};    // Intermediate search direction
  vector<double> mult{}; // Temporary vector for scalar multiplication
};

// Represents boundary conditions.
// Stores the type, function identifier, and node numbers on the boundary.
struct bc
{
  int type{};          // Type of boundary condition (e.g., 1 for Dirichlet)
  int function{};      // Identifier for the boundary condition function
  vector<int> ndnum{}; // Global node numbers on the boundary
};

/*============================== I/O  FUNCTIONS =============================*/

// Reads mesh nodes, elements, properties, and function number from files.
void input(vector<nd> &mesh, vector<el> &elList, double &sigma, int &fnum);
// Reads coarse time steps from a file.
void input_timemesh(TMesh &time);
// Refines the time mesh by splitting intervals based on input data.
void input_split_timemesh(TMesh &time);
// Reads boundary condition data from a file.
void input_boundary(vector<bc> &cond);
// Prints the solution to the terminal for each time step and node.
void print_solution(const TMesh &timemesh, const vector<nd> &mesh);
// Outputs the solution to a file, including numerical and analytical solutions
// and errors.
void out_solution(TMesh &timemesh, vector<nd> &mesh, int flag);

/*========================= LOCAL MATRICES FUNCTIONS =======================*/

// Computes the determinant of the element's coordinate matrix.
double detD(el &e);
// Computes the factorial of a number.
long fact(int a);
// Integrates the product of shape functions L1*L2*L3.
double intL(int nv[], double det);
// Computes an element of the mass matrix M_ij.
double Mij(int i, int j, el e, double det);
// Computes the local mass matrix M.
void getM(vector<vector<double>> &M, double sigma, el &e);
// Computes the local stiffness matrix G.
void getG(vector<vector<double>> &G, el e);
// Computes the local load vector b.
void getb(vector<double> &b, el e, double t, double gamma, int flag);

/*======================== SPARSE MATRICES FUNCTIONS =======================*/
// Constructs the sparse matrix structure (portrait).
void portrait(vector<nd> &mesh, vector<el> &elList, SLAE &slae);
// Adds a local matrix to the global sparse matrix.
void add_mx(smx &A, el elem, vector<vector<double>> &localMatrix);
// Adds a local vector to the global vector.
void add_vec(vector<double> &b, el &elem, vector<double> &bLocal);
// Adds a single element to the sparse matrix at position (i, j).
void add_el(smx &A, int i, int j, double elem);
// Multiplies a local matrix by a global vector.
void mult_mx_vec(vector<vector<double>> &matrix, vector<double> &vec,
                 vector<double> &result, el &elem);
// Multiplies a matrix by a scalar.
void mult_mx_num(vector<vector<double>> &matrix, double coef,
                 vector<vector<double>> &resultMatrix);
// Multiplies a vector by a scalar in-place.
void mult_vec_num(vector<double> &vector, double coef);

/*============================== LOS FUNCTIONS =============================*/

// Performs LU factorization of the sparse matrix.
void calcLU(SLAE &slae, SLAE &LU);
// Solves the system using the Conjugate Gradient method with LU
// preconditioning.
void losLU(SLAE &slae, SLAE &LU, LOS &v, int maxIter, double eps);
// Solves Ly = b using forward substitution.
void calcY(SLAE &LU, vector<double> &b, vector<double> &y);
// Solves Ux = y using backward substitution.
void calcX(SLAE &LU, vector<double> &y, vector<double> &x);
// Multiplies a sparse matrix by a vector.
void mult_smx_vec(smx &A, vector<double> &x, vector<double> &F);
// Computes the residual (discrepancy) of the system Ax = b.
void calc_discrepancy(SLAE &slae, LOS &v, vector<double> &x, double &normb);
// Clears the sparse matrix and vectors in the SLAE structure.
void clearSLAE(SLAE &slae);

/*============================= SOLVER FUNCTIONS =============================*/

// Main solver function for the time-dependent problem.
void Solver(vector<nd> &mesh, vector<el> &elList, TMesh &timemesh, SLAE &slae,
            vector<bc> &conds, int flag, double sigma);
// Assembles the global system matrix and vector for a given time step.
void get_global(vector<nd> &mesh, vector<el> elList, TMesh &timemesh,
                SLAE &slae, vector<bc> &cond, int ti, int flag, double sigma);
// Initializes solution vectors for the first two time steps.
void initq0q1(TMesh &timemesh, vector<nd> &mesh, int flag);
// Applies first-type (Dirichlet) boundary conditions.
void bc1(SLAE &slae, vector<bc> &cond, vector<nd> mesh, double tValue);
// Clears a row in the sparse matrix, leaving only the diagonal element.
void mx_clearline(smx &A, int i);

/*============================== F/U FUNCTIONS =============================*/

// Defines the source term f(n, t) based on the flag.
double f(nd &n, double t, int flag);
// Defines the analytical solution u(n, t) based on the flag.
double u(nd &n, double t, int flag);