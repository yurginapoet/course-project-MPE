#include "header.h"

// Initializes the solution vectors qti_2 and qti_1 for the first two time
// steps. Uses the analytical solution u to set initial conditions at t =
// tNew[0] and t = tNew[1].
void initq0q1(TMesh &timemesh, vector<nd> &mesh, int flag)
{
  auto &qti = timemesh.qti;     // Current solution vector
  auto &qti_1 = timemesh.qti_1; // Solution at previous time step
  auto &qti_2 = timemesh.qti_2; // Solution at two time steps back

  int mxsize = mesh.size(); // Number of nodes

  qti.resize(mxsize);   // Resize current solution
  qti_1.resize(mxsize); // Resize previous solution
  qti_2.resize(mxsize); // Resize solution at t-2

  // Set initial conditions using the analytical solution u
  for (int i = 0; i < mxsize; i++)
  {
    qti_2[i] = u(mesh[i], timemesh.tNew[0], flag); // Solution at t = tNew[0]
    qti_1[i] = u(mesh[i], timemesh.tNew[1], flag); // Solution at t = tNew[1]
  }
}

// Assembles the global system matrix and right-hand side vector for a given
// time step. Incorporates mass matrix M, stiffness matrix G, and load vector b,
// scaled by time step coefficients.
void get_global(vector<nd> &mesh, vector<el> elList, TMesh &timemesh,
                SLAE &slae, vector<bc> &cond, int ti, int flag, double sigma)
{
  auto &A = slae.A; // Global sparse matrix
  auto &b = slae.b; // Global right-hand side vector

  double t = timemesh.tNew[ti]; // Current time
  double deltaT =
      timemesh.tNew[ti] - timemesh.tNew[ti - 2]; // Time step from t-2 to t
  double deltaT0 =
      timemesh.tNew[ti] - timemesh.tNew[ti - 1]; // Time step from t-1 to t
  double deltaT1 = timemesh.tNew[ti - 1] -
                   timemesh.tNew[ti - 2]; // Time step from t-2 to t-1
  int slaesize = mesh.size();             // Number of nodes

  vector<vector<double>> M(3), G(3),
      tempMatrix(3); // Local mass, stiffness, and temporary matrices
  vector<double> locb(3, 0),
      tempVector(3, 0); // Local load and temporary vectors

  // Resize local matrices
  for (int i = 0; i < 3; i++)
  {
    M[i].resize(3);
    G[i].resize(3);
    tempMatrix[i].resize(3);
  }

  // Loop over all elements to assemble the global system
  for (auto &elem : elList)
  {
    // Calculate local matrices and vectors
    getG(G, elem);                    // Compute stiffness matrix
    getM(M, sigma, elem);             // Compute mass matrix
    getb(locb, elem, t, sigma, flag); // Compute local load vector

    // Scale mass matrix and add to global matrix
    mult_mx_num(M, (deltaT + deltaT0) / (deltaT * deltaT0),
                tempMatrix);     // Scale M
    add_mx(A, elem, tempMatrix); // Add scaled M to global matrix
    add_mx(A, elem, G);          // Add G to global matrix

    // Add local load vector to global vector
    add_vec(b, elem, locb);

    // Add contributions from previous time steps
    fill(tempVector.begin(), tempVector.end(), 0.0);
    mult_vec_num(tempVector, -deltaT0 / (deltaT * deltaT1)); // Scale for qti_2
    mult_mx_vec(M, timemesh.qti_2, tempVector, elem);        // M * qti_2
    add_vec(b, elem, tempVector);                            // Add to global b

    fill(tempVector.begin(), tempVector.end(), 0.0);
    mult_mx_vec(M, timemesh.qti_1, tempVector, elem);       // M * qti_1
    mult_vec_num(tempVector, deltaT / (deltaT1 * deltaT0)); // Scale for qti_1
    add_vec(b, elem, tempVector);                           // Add to global b
  }
  // Apply first-type boundary conditions
  bc1(slae, cond, mesh, t);
}

// Applies first-type (Dirichlet) boundary conditions to the system.
// Sets the corresponding rows in the matrix to represent fixed values.
void bc1(SLAE &slae, vector<bc> &cond, vector<nd> mesh, double tValue)
{
  auto &A = slae.A; // Sparse matrix
  auto &b = slae.b; // Right-hand side vector

  // Loop over boundary conditions
  for (auto &condition : cond)
    if (condition.type == 1) // First-type (Dirichlet) condition
    {
      auto &ndnum = condition.ndnum; // Node indices for boundary
      int flag = condition.function; // Function identifier

      const int numVertex = ndnum.size(); // Number of boundary nodes

      for (int i = 0; i < numVertex; i++)
      {
        int globalNum = ndnum[i];   // Global node number
        mx_clearline(A, globalNum); // Clear row, leaving only diagonal
        b[globalNum] = u(mesh[globalNum], tValue, flag); // Set boundary value
      }
    }
}

// Clears a row in the sparse matrix A, setting the diagonal to 1 and
// off-diagonal elements to 0. Used for applying Dirichlet boundary conditions.
void mx_clearline(smx &A, int i)
{
  auto &ig = A.ig, &jg = A.jg; // Row pointers and column indices
  auto &ggl = A.ggl, &ggu = A.ggu, &di = A.di; // Matrix components
  const int size = ig[di.size()]; // Total number of off-diagonal elements

  int i0 = ig[i], i1 = ig[i + 1]; // Range for row i

  // Clear lower triangular elements
  for (i0; i0 < i1; i0++)
    ggl[i0] = 0.;

  // Clear upper triangular elements where column index matches i
  int j0 = i1;
  for (j0; j0 < size; j0++)
    if (jg[j0] == i)
      ggu[j0] = 0.;

  di[i] = 1.; // Set diagonal to 1
}

// Main solver function for the time-dependent problem.
// Iterates over time steps, assembles and solves the system, and updates
// solutions.
void Solver(vector<nd> &mesh, vector<el> &elList, TMesh &timemesh, SLAE &slae,
            vector<bc> &conds, int flag, double sigma)
{
  // Initialize solutions for the first two time steps
  initq0q1(timemesh, mesh, flag);

  // Store initial solutions
  timemesh.q[0] = timemesh.qti_2; // t = tNew[0]
  timemesh.q[1] = timemesh.qti_1; // t = tNew[1]

  int tSize = timemesh.tNew.size();  // Number of time steps
  for (int ti = 2; ti < tSize; ++ti) // Start from third time step
  {
    double deltaT =
        timemesh.tNew[ti] - timemesh.tNew[ti - 2]; // Time step from t-2 to t
    double deltaT0 =
        timemesh.tNew[ti] - timemesh.tNew[ti - 1]; // Time step from t-1 to t
    double deltaT1 = timemesh.tNew[ti - 1] -
                     timemesh.tNew[ti - 2]; // Time step from t-2 to t-1

    SLAE LU{};                    // LU factorization structure
    LOS vectors{};                // LOS vectors for iterative solver
    auto &A = slae.A;             // Global sparse matrix
    auto &b = slae.b;             // Right-hand side vector
    double t = timemesh.tNew[ti]; // Current time
    int slaesize = mesh.size();   // Number of nodes

    // Initialize local matrices and vectors
    vector<vector<double>> M(3, vector<double>(3, 0.0)); // Mass matrix
    vector<vector<double>> G(3, vector<double>(3, 0.0)); // Stiffness matrix
    vector<vector<double>> tempMatrix(
        3, vector<double>(3, 0.0)); // Temporary matrix
    vector<double> locb(3, 0.0),
        tempVector(3, 0.0); // Local load and temporary vectors

    // Assemble the global system
    for (size_t e = 0; e < elList.size(); ++e)
    {
      auto &elem = elList[e];           // Current element
      getM(M, sigma, elem);             // Compute mass matrix
      getG(G, elem);                    // Compute stiffness matrix
      getb(locb, elem, t, sigma, flag); // Compute local load vector

      double coef =
          (deltaT + deltaT0) / (deltaT * deltaT0); // Time step coefficient
      mult_mx_num(M, coef, tempMatrix);            // Scale mass matrix
      add_mx(A, elem, tempMatrix);                 // Add to global matrix
      add_mx(A, elem, G);                          // Add stiffness matrix
      add_vec(b, elem, locb);                      // Add local load vector

      // Add contributions from previous time steps
      fill(tempVector.begin(), tempVector.end(), 0.0);
      mult_mx_vec(M, timemesh.qti_2, tempVector, elem);        // M * qti_2
      mult_vec_num(tempVector, -deltaT0 / (deltaT * deltaT1)); // Scale
      add_vec(b, elem, tempVector); // Add to global b

      fill(tempVector.begin(), tempVector.end(), 0.0);
      mult_mx_vec(M, timemesh.qti_1, tempVector, elem);       // M * qti_1
      mult_vec_num(tempVector, deltaT / (deltaT1 * deltaT0)); // Scale
      add_vec(b, elem, tempVector);                           // Add to global b
    }

    // Apply boundary conditions
    bc1(slae, conds, mesh, t);
    // Perform LU factorization
    calcLU(slae, LU);
    // Solve the system using LOS with LU preconditioning
    losLU(slae, LU, vectors, 10000, 1e-30);

    // Store the solution and update time step vectors
    timemesh.qti = slae.q; // Current solution
    timemesh.save_ti(ti);  // Save solution for current time step
    timemesh.Swap();       // Update qti_2 and qti_1 for next iteration
    clearSLAE(slae);       // Clear the system for the next time step
  }
}