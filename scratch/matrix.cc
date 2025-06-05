#include "header.h"

/*======================= MATRIX & VECTOR OPERATIONS =========================*/

// Adds a local matrix to the global sparse matrix A.
// Iterates over the 3x3 local matrix and adds each element to the appropriate
// position in the sparse matrix structure (diagonal, lower triangular, or upper
// triangular).
void add_mx(smx &A, el elem, vector<vector<double>> &localMatrix)
{
  const int sizeLocal = 3; // Local matrix is 3x3 (for triangular elements)
  auto locV = elem.nds;    // Nodes of the current element

  // Loop through the local matrix
  for (int i = 0; i < sizeLocal; i++)
    for (int j = 0; j < sizeLocal; j++)
    {
      double elem = localMatrix[i][j]; // Value from the local matrix
      // Add the element to the global sparse matrix
      add_el(A, locV[i].gl_num, locV[j].gl_num, elem);
    }
}

// Adds a local vector to the global vector b.
// Accumulates contributions from the local vector to the global vector at
// positions corresponding to the global node numbers.
void add_vec(vector<double> &b, el &elem, vector<double> &bLocal)
{
  const int sizeLocal =
      bLocal.size();    // Size of the local vector (3 for triangular elements)
  auto locV = elem.nds; // Nodes of the current element

  // Add local vector contributions to the global vector
  for (int i = 0; i < sizeLocal; i++)
    b[locV[i].gl_num] += bLocal[i];
}

// Adds a single element to the sparse matrix A at position (i, j).
// Handles diagonal (i == j), lower triangular (i > j), and upper triangular (i
// < j) cases. Uses binary search to locate the correct position in the sparse
// matrix structure.
void add_el(smx &A, int i, int j, double elem)
{
  auto &ig = A.ig, &jg = A.jg; // Row pointers and column indices
  auto &ggl = A.ggl, &ggu = A.ggu,
       &di = A.di; // Lower, upper, and diagonal elements

  if (i == j) // Diagonal element
    di[i] += elem;
  else if (i > j) // Lower triangular element
  {
    int beg = ig[i], end = ig[i + 1] - 1; // Range for binary search
    while (jg[beg] != j)
    {
      int ind = (beg + end) / 2; // Middle index
      if (jg[ind] < j)
        beg = ind + 1;
      else
        end = ind;
    }
    ggl[beg] += elem; // Add to lower triangular part
  }
  else // Upper triangular element
  {
    int beg = ig[j], end = ig[j + 1] - 1; // Range for binary search
    while (jg[beg] != i)
    {
      int ind = (beg + end) / 2; // Middle index
      if (jg[ind] < i)
        beg = ind + 1;
      else
        end = ind;
    }
    ggu[beg] += elem; // Add to upper triangular part
  }
}

// Multiplies a local matrix by a global vector, producing a local result
// vector. The global vector values are accessed using global node numbers of
// the element.
void mult_mx_vec(vector<vector<double>> &matrix, vector<double> &vec,
                 vector<double> &result, el &elem)
{
  const int sizeMatrix = matrix.size(); // Size of local matrix (3x3)
  auto locV = elem.nds;                 // Nodes of the element

  for (int i = 0; i < sizeMatrix; i++)
  {
    double sum = 0; // Accumulate result for row i
    for (int j = 0; j < sizeMatrix; j++)
      sum += matrix[i][j] * vec[locV[j].gl_num]; // Multiply and accumulate
    result[i] = sum;                             // Store result
  }
}

// Multiplies a matrix by a scalar coefficient.
// The result is stored in a provided result matrix.
void mult_mx_num(vector<vector<double>> &matrix, double coef,
                 vector<vector<double>> &resultMatrix)
{
  const int sizeMatrix = matrix.size(); // Size of matrix (3x3)

  for (int i = 0; i < sizeMatrix; i++)
    for (int j = 0; j < sizeMatrix; j++)
      resultMatrix[i][j] = coef * matrix[i][j]; // Scale each element
}

// Multiplies a vector by a scalar coefficient in-place.
void mult_vec_num(vector<double> &vector, double coef)
{
  const int size = vector.size(); // Vector size

  for (int i = 0; i < size; i++)
    vector[i] *= coef; // Scale each element
}

// Computes the sum of two vectors a and b, storing the result in res.
void sum_vec(vector<double> &a, vector<double> &b, vector<double> &res)
{
  const int n = a.size(); // Vector size

  for (int i = 0; i < n; i++)
    res[i] = a[i] + b[i]; // Element-wise sum
}

// Computes the dot product of two vectors a and b.
double dotproduct(vector<double> &a, vector<double> &b)
{
  int n = a.size(); // Vector size
  double res = 0;
  for (int i = 0; i < n; i++)
    res += a[i] * b[i]; // Accumulate product
  return res;
}

// Multiplies vector a by a scalar coef, storing the result in res.
void mult_vec_num(vector<double> &a, double coef, vector<double> &res)
{
  const int n = a.size(); // Vector size

  for (int i = 0; i < n; i++)
    res[i] = a[i] * coef; // Scale each element
}

// Multiplies a sparse matrix A by a vector x, storing the result in F.
// Handles diagonal, lower, and upper triangular contributions.
void mult_smx_vec(smx &A, vector<double> &x, vector<double> &F)
{
  auto &ig = A.ig, &jg = A.jg; // Row pointers and column indices
  auto &di = A.di, &ggl = A.ggl,
       &ggu = A.ggu;           // Diagonal, lower, upper elements
  const int sizeA = di.size(); // Size of the matrix

  for (int i = 0; i < sizeA; i++)
  {
    F[i] = di[i] * x[i];            // Diagonal contribution
    int i0 = ig[i], i1 = ig[i + 1]; // Range of non-zero elements in row i

    for (i0; i0 < i1; i0++) // Process lower and upper triangular parts
    {
      int j = jg[i0];         // Column index
      F[i] += ggl[i0] * x[j]; // Lower triangular contribution
      F[j] += ggu[i0] * x[i]; // Upper triangular contribution
    }
  }
}

// Computes the discrepancy (residual) of the system Ax = b.
// Stores the result in v.r1 and computes the norm of b.
void calc_discrepancy(SLAE &slae, LOS &v, vector<double> &x, double &normb)
{
  auto &ig = slae.A.ig, &jg = slae.A.jg; // Row pointers and column indices
  auto &ggl = slae.A.ggl, &ggu = slae.A.ggu, &di = slae.A.di, &b = slae.b,
       &r1 = v.r1;                // Matrix components and residual vector
  const int sizeSlae = di.size(); // Size of the system

  for (int i = 0; i < sizeSlae; i++)
  {
    normb += b[i] * b[i];           // Compute norm of b
    r1[i] = b[i] - di[i] * x[i];    // Diagonal contribution to residual
    int i0 = ig[i], i1 = ig[i + 1]; // Range of non-zero elements
    for (i0; i0 < i1; i0++)
    {
      int j = jg[i0];          // Column index
      r1[i] -= ggl[i0] * x[j]; // Lower triangular contribution
      r1[j] -= ggu[i0] * x[i]; // Upper triangular contribution
    }
  }
}

// Clears the sparse matrix A and vectors b and q in the SLAE structure.
// Resets all elements to zero.
void clearSLAE(SLAE &slae)
{
  auto &A = slae.A; // Sparse matrix
  auto &b = slae.b; // Right-hand side vector
  auto &q = slae.q; // Solution vector

  const int sizeSLAE = A.di.size(); // Size of the system
  const int countNonZeroElems =
      A.ig[sizeSLAE]; // Number of non-zero off-diagonal elements

  for (int i = 0; i < sizeSLAE; i++)
  {
    A.di[i] = 0; // Clear diagonal
    q[i] = 0;    // Clear solution
    b[i] = 0;    // Clear right-hand side
  }

  for (int i = 0; i < countNonZeroElems; i++)
  {
    A.ggl[i] = 0; // Clear lower triangular
    A.ggu[i] = 0; // Clear upper triangular
  }
}

// Constructs the portrait (structure) of the sparse matrix.
// Determines the non-zero pattern of the matrix based on element connectivity.
void portrait(vector<nd> &mesh, vector<el> &elList, SLAE &slae)
{
  auto &nodeCoord = mesh;  // Mesh nodes
  auto &elements = elList; // Element list

  auto &ig = slae.A.ig, &jg = slae.A.jg; // Row pointers and column indices
  auto &b = slae.b;                      // Right-hand side vector

  vector<set<int>> rowCount{}; // Temporary structure to store non-zero column
                               // indices per row

  const int sizeSlae = nodeCoord.size(); // Number of nodes (matrix size)

  slae.A.di.resize(sizeSlae); // Resize diagonal vector
  slae.b.resize(sizeSlae);    // Resize right-hand side vector
  slae.q.resize(sizeSlae);    // Resize solution vector
  ig.resize(sizeSlae + 1);    // Resize row pointers
  rowCount.resize(sizeSlae);  // Resize temporary structure

  // Build the non-zero structure by checking element connectivity
  for (auto &elem : elements)
  {
    for (auto &i : elem.nds)
    {
      for (auto &j : elem.nds)
        if (j.gl_num < i.gl_num) // Store only lower triangular indices
          rowCount[i.gl_num].insert(j.gl_num);
    }
  }

  ig[0] = 0; // Initialize first row pointer
  for (int i = 0; i < sizeSlae; i++)
  {
    ig[i + 1] = ig[i] + rowCount[i].size(); // Set row pointers
    for (auto j : rowCount[i])
      jg.push_back(j); // Add column indices
  }

  slae.A.ggl.resize(ig[sizeSlae]); // Resize lower triangular elements
  slae.A.ggu.resize(ig[sizeSlae]); // Resize upper triangular elements
}

/*=================================== LOS ====================================*/

// Solves the system Ax = b using the Conjugate Gradient method with LU
// preconditioning. Iteratively refines the solution until convergence or max
// iterations is reached.
void losLU(SLAE &slae, SLAE &LU, LOS &v, int maxIter, double eps)
{
  auto &r1 = v.r1, &z1 = v.z1, &p1 = v.p1, &mult = v.mult, &rk = v.rk,
       &Ar = v.Ar, &p = v.p;      // LOS vectors
  auto &q = slae.q;               // Solution vector
  const int n = slae.A.di.size(); // System size
  double normb = 0;               // Norm of the right-hand side

  // Resize working vectors
  p.resize(n);
  r1.resize(n);
  z1.resize(n);
  p1.resize(n);
  mult.resize(n);
  rk.resize(n);
  Ar.resize(n);

  // Compute initial residual r = b - Ax
  calc_discrepancy(slae, v, q, normb);
  if (abs(normb) < 1e-10) // Check for invalid norm
  {
    cerr << "Error: normb is zero or too small: " << normb << endl;
    exit(1);
  }

  // Solve LUz = r for initial z
  calcY(LU, r1, r1);                            // Forward substitution
  calcX(LU, r1, z1);                            // Backward substitution
  mult_smx_vec(slae.A, z1, p1);                 // Compute p1 = Az1
  calcY(LU, p1, p1);                            // Solve LUp1 = Az1
  double scalarr = dotproduct(r1, r1);          // Compute ||r||^2
  double discrepancy = sqrt(scalarr / normb);   // Initial discrepancy
  if (isnan(discrepancy) || isinf(discrepancy)) // Check for numerical issues
  {
    cerr << "Error: Initial discrepancy is NaN or Inf" << endl;
    exit(1);
  }

  // Iterative refinement using Conjugate Gradient
  for (int k = 1; k < maxIter && discrepancy > eps; k++)
  {
    double scalarp = dotproduct(p1, p1); // Compute ||p1||^2
    if (abs(scalarp) < 1e-35)            // Check for numerical stability
    {
      cerr << "Error: scalarp is zero or too small at iteration " << k
           << "it equals:" << scalarp << endl;
      exit(1);
    }
    double alpha = dotproduct(p1, r1) / scalarp; // Compute step size
    mult_vec_num(z1, alpha, mult);               // mult = alpha * z1
    sum_vec(q, v.mult, q);          // Update solution: q = q + alpha * z1
    mult_vec_num(p1, -alpha, mult); // mult = -alpha * p1
    sum_vec(r1, mult, r1);          // Update residual: r1 = r1 - alpha * p1

    calcX(LU, r1, rk);                           // Solve LUx = r1 for rk
    mult_smx_vec(slae.A, rk, Ar);                // Ar = A * rk
    calcY(LU, Ar, p);                            // Solve LUp = Ar
    double betta = -dotproduct(p1, p) / scalarp; // Compute beta
    mult_vec_num(z1, betta, mult);               // mult = beta * z1
    sum_vec(rk, mult, z1);         // Update z1: z1 = rk + beta * z1
    mult_vec_num(p1, betta, mult); // mult = beta * p1
    sum_vec(p, mult, p1);          // Update p1: p1 = p + beta * p1
    discrepancy = sqrt(dotproduct(r1, r1) / normb); // Update discrepancy
    if (isnan(discrepancy) || isinf(discrepancy)) // Check for numerical issues
    {
      cerr << "Error: Discrepancy is NaN or Inf at iteration " << k << endl;
      exit(1);
    }
  }
  // Final discrepancy check
  normb = 0;
  calc_discrepancy(slae, v, q, normb);
  discrepancy = sqrt(dotproduct(r1, r1) / normb);
}

// Performs LU factorization of the sparse matrix A, storing L and U in LU.
void calcLU(SLAE &slae, SLAE &LU)
{
  auto &ig = slae.A.ig, &jg = slae.A.jg; // Row pointers and column indices
  auto &ggl = slae.A.ggl, &ggu = slae.A.ggu,
       &di = slae.A.di;                              // Matrix components
  auto &L = LU.A.ggl, &U = LU.A.ggu, &diL = LU.A.di; // LU components
  LU.b = slae.b;                                     // Copy right-hand side
  LU.A.ig = ig;                                      // Copy row pointers
  LU.A.jg = jg;                                      // Copy column indices

  const int sizeSlae = di.size(), sizeTriangles = ig[sizeSlae]; // Sizes

  diL.resize(sizeSlae, 0.0);    // Initialize diagonal of L
  L.resize(sizeTriangles, 0.0); // Initialize lower triangular part
  U.resize(sizeTriangles, 0.0); // Initialize upper triangular part

  for (int i = 0; i < sizeSlae; i++)
  {
    double sumDi = 0.0;             // Accumulate diagonal contribution
    int i0 = ig[i], i1 = ig[i + 1]; // Range of non-zero elements in row i

    // Compute L and U elements for row i
    for (int k = i0; k < i1; k++)
    {
      int j = jg[k];                 // Column index
      double suml = 0.0, sumu = 0.0; // Accumulators for L and U

      // Compute contributions from previous rows
      for (int m = ig[j]; m < ig[j + 1]; m++)
      {
        for (int n = ig[i]; n < ig[i + 1]; n++)
        {
          if (jg[m] == jg[n] && jg[m] < j) // Match column indices
          {
            suml += L[n] * U[m]; // Contribution to L
            sumu += L[m] * U[n]; // Contribution to U
          }
        }
      }

      L[k] = ggl[k] - suml;    // Compute L element
      if (abs(diL[j]) < 1e-10) // Check for near-zero diagonal
      {
        cerr << "Error: Near-zero diagonal element diL[" << j
             << "] = " << diL[j] << endl;
        exit(1);
      }
      U[k] = (ggu[k] - sumu) / diL[j]; // Compute U element
      sumDi += L[k] * U[k];            // Accumulate for diagonal
    }

    // Compute diagonal element of L
    diL[i] = di[i] - sumDi;
    if (abs(diL[i]) < 1e-10) // Check for near-zero diagonal
    {
      cerr << "Error: Near-zero diagonal element diL[" << i << "] = " << diL[i]
           << endl;
      exit(1);
    }
  }
}

// Solves Ly = b using forward substitution (L is lower triangular).
void calcY(SLAE &LU, vector<double> &b, vector<double> &y)
{
  auto &ig = LU.A.ig, &jg = LU.A.jg; // Row pointers and column indices
  auto &di = LU.A.di, &L = LU.A.ggl; // Diagonal and lower triangular parts
  const int sizeSlae = di.size();    // System size

  for (int i = 0; i < sizeSlae; i++)
  {
    double sum = 0;                 // Accumulate contributions
    int i0 = ig[i], i1 = ig[i + 1]; // Range of non-zero elements

    for (i0; i0 < i1; i0++)
    {
      int j = jg[i0];      // Column index
      sum += L[i0] * y[j]; // Contribution from previous y values
    }

    y[i] = (b[i] - sum) / di[i]; // Compute y[i]
  }
}

// Solves Ux = y using backward substitution (U is upper triangular).
void calcX(SLAE &LU, vector<double> &y, vector<double> &x)
{
  auto &ig = LU.A.ig, &jg = LU.A.jg;   // Row pointers and column indices
  auto &U = LU.A.ggu;                  // Upper triangular part
  const int sizeSlae = LU.A.di.size(); // System size
  vector<double> v = y;                // Temporary vector

  for (int i = sizeSlae - 1; i >= 0; i--)
  {
    x[i] = v[i];                    // Set x[i]
    int i0 = ig[i], i1 = ig[i + 1]; // Range of non-zero elements

    for (i0; i0 < i1; i0++)
    {
      int j = jg[i0];       // Column index
      v[j] -= x[i] * U[i0]; // Update v[j] for next iterations
    }
  }
}