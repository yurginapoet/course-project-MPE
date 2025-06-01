#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace std;

// Structure for mesh node
struct nd
{
  int gl_num;
  double r, z;
  nd(int _num, double _r, double _z) : gl_num(_num), r(_r), z(_z) {}
  nd() = default;
};

// Structure for element
struct el
{
  double lambda;
  nd nds[3];
};

// Input data
void input(vector<nd> &mesh, vector<el> &elList, double &gamma, int &fnum,
           vector<double> &time)
{
  int n, node;
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

  // Read time steps
  in.open("../data/time.txt");
  in >> n;
  time.resize(n);
  for (int i = 0; i < n; i++)
  {
    in >> time[i];
  }
  in.close();
}

// Output results
void output(const vector<double> &q, const vector<nd> &mesh, double t)
{
  ofstream out("data/q.txt", ios::app);
  out << scientific << setprecision(6);
  out << "Time: " << t << "\n";
  for (int i = 0; i < q.size(); ++i)
  {
    out << "node " << mesh[i].gl_num << " (r = " << mesh[i].r
        << ", z = " << mesh[i].z << "): q = " << q[i] << "\n";
  }
  out << "\n";
  out.close();
}

// Compute determinant
double detD(const el &e)
{
  return abs((e.nds[1].r - e.nds[0].r) * (e.nds[2].z - e.nds[0].z) -
             (e.nds[2].r - e.nds[0].r) * (e.nds[1].z - e.nds[0].z));
}

// Factorial
long fact(int a)
{
  long f = 1;
  for (int i = 1; i <= a; i++)
    f *= i;
  return f;
}

// Integrate L1*L2*L3
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

// Compute mass matrix element M_ij
double Mij(int i, int j, const el &e, double det)
{
  int nv[3];
  double sum = 0;
  for (int m = 0; m < 3; m++)
  {
    nv[0] = nv[1] = nv[2] = 0;
    nv[i]++;
    nv[j]++;
    nv[m]++;
    sum += e.nds[m].r * intL(nv, det);
  }
  return sum;
}

// Get local mass matrix
void getM(vector<vector<double>> &M, double gamma, const el &e)
{
  double det = detD(e);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      M[i][j] = gamma * Mij(i, j, e, det);
}

// Get local stiffness matrix
void getG(vector<vector<double>> &G, const el &e)
{
  vector<double> a1(3), a2(3);
  a1[0] = (e.nds[1].z - e.nds[2].z);
  a1[1] = (e.nds[2].z - e.nds[0].z);
  a1[2] = (e.nds[0].z - e.nds[1].z);
  a2[0] = (e.nds[2].r - e.nds[1].r);
  a2[1] = (e.nds[0].r - e.nds[2].r);
  a2[2] = (e.nds[1].r - e.nds[0].r);

  double det =
      (detD(e) * e.lambda) * (e.nds[0].r + e.nds[1].r + e.nds[2].r) / 12;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      G[i][j] = det * (a1[i] * a1[j] + a2[i] * a2[j]);
}

// Get local load vector
void getb(vector<double> &b, const el &e, double sigma, int s, double t,
          const vector<double> &q1, const vector<double> &q2, double dt)
{
  double det = detD(e);
  vector<double> p_loc(3);
  for (int i = 0; i < 3; i++)
  {
    p_loc[i] = f(s, e.nds[i], e); // Assuming f is defined as in your code
  }

  vector<vector<double>> M(3, vector<double>(3));
  getM(M, sigma, e);

  vector<double> q1_loc(3), q2_loc(3);
  for (int i = 0; i < 3; i++)
  {
    q1_loc[i] = q1[e.nds[i].gl_num];
    q2_loc[i] = q2[e.nds[i].gl_num];
  }

  for (int i = 0; i < 3; i++)
  {
    b[i] = 0;
    for (int j = 0; j < 3; j++)
    {
      b[i] += M[i][j] * (p_loc[j] + (4.0 * q1_loc[j] - q2_loc[j]) / (2.0 * dt));
    }
  }
}

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

// Build sparse matrix portrait
void buildPortrait(vector<vector<int>> &list, vector<int> &ig, vector<int> &jg,
                   const vector<el> &elList)
{
  int Nuz = list.size();
  list[0].push_back(0);
  for (const auto &e : elList)
  {
    for (int i = 0; i < 3; i++)
    {
      for (int j = i + 1; j < 3; j++)
      {
        int insertPos = e.nds[j].gl_num;
        int element = e.nds[i].gl_num;
        bool isIn = false;
        for (int k : list[insertPos])
        {
          if (k == element)
          {
            isIn = true;
            break;
          }
        }
        if (!isIn)
          list[insertPos].push_back(element);
      }
    }
  }

  for (int i = 0; i < Nuz; i++)
  {
    sort(list[i].begin(), list[i].end());
  }

  ig[0] = 0;
  ig[1] = 0;
  for (int i = 1; i < Nuz; i++)
  {
    ig[i + 1] = ig[i] + list[i].size();
  }

  jg.resize(ig.back());
  for (int i = 1, j = 0; i < Nuz; i++)
  {
    for (int k = 0; k < list[i].size(); k++)
    {
      jg[j++] = list[i][k];
    }
  }
}

// Add element to sparse matrix
void addElementToGlobal(int i, int j, double elem, vector<double> &di,
                        vector<double> &gg, const vector<int> &ig,
                        const vector<int> &jg)
{
  if (i == j)
  {
    di[i] += elem;
    return;
  }
  for (int ind = ig[i]; ind < ig[i + 1]; ind++)
  {
    if (jg[ind] == j)
    {
      gg[ind] += elem;
      return;
    }
  }
}

// Apply first boundary conditions
void applyBC1(vector<double> &di, vector<double> &gg, vector<double> &F,
              const vector<int> &ig, const vector<int> &jg,
              const vector<nd> &mesh)
{
  ifstream in("../data/c1.txt");
  int p;
  in >> p;
  vector<int> bc1nodes(mesh.size(), -1);
  vector<pair<int, double>> bc1(p * 2);
  for (int i = 0; i < p; i++)
  {
    int v1, v2, s;
    in >> v1 >> v2 >> s;
    bc1[2 * i] = {v1, u(s, mesh[v1])};
    bc1[2 * i + 1] = {v2, u(s, mesh[v2])};
    bc1nodes[v1] = 2 * i;
    bc1nodes[v2] = 2 * i + 1;
  }
  in.close();

  for (int i = 0; i < mesh.size(); i++)
  {
    if (bc1nodes[i] != -1)
    {
      di[i] = 1.0;
      F[i] = bc1[bc1nodes[i]].second;
      for (int j = ig[i]; j < ig[i + 1]; j++)
      {
        int k = jg[j];
        if (bc1nodes[k] == -1)
        {
          F[k] -= gg[j] * F[i];
        }
        gg[j] = 0.0;
      }
    }
    else
    {
      for (int j = ig[i]; j < ig[i + 1]; j++)
      {
        int k = jg[j];
        if (bc1nodes[k] != -1)
        {
          F[k] -= gg[j] * F[i];
          gg[j] = 0.0;
        }
      }
    }
  }
}

// Matrix-vector multiplication for sparse matrix
void mult(const vector<double> &di, const vector<double> &gg,
          const vector<int> &ig, const vector<int> &jg, const vector<double> &x,
          vector<double> &y)
{
  for (int i = 0; i < y.size(); i++)
  {
    y[i] = di[i] * x[i];
    for (int k = ig[i]; k < ig[i + 1]; k++)
    {
      int j = jg[k];
      y[i] += gg[k] * x[j];
      y[j] += gg[k] * x[i];
    }
  }
}

// Euclidean norm
double euclideanNorm(const vector<double> &x)
{
  double scalar = 0;
  for (double val : x)
  {
    scalar += val * val;
  }
  return sqrt(scalar);
}

// Conjugate gradient solver
void conjugateGradient(const vector<double> &di, const vector<double> &gg,
                       const vector<int> &ig, const vector<int> &jg,
                       const vector<double> &F, vector<double> &q, double eps,
                       int maxIter)
{
  int Nuz = q.size();
  vector<double> um(Nuz, 0), r(Nuz, 0), z(Nuz, 0);
  mult(di, gg, ig, jg, q, um);
  for (int i = 0; i < Nuz; i++)
  {
    r[i] = F[i] - um[i];
  }
  z = r;
  double bnorm = euclideanNorm(F);
  double residual = euclideanNorm(r) / bnorm;
  if (residual <= eps)
    return;

  double scal1 = 0, scal2 = 0, scal3 = 0, alfa = 0, beta = 0;
  mult(z, um);
  for (int i = 0; i < Nuz; i++)
  {
    scal1 += r[i] * r[i];
    scal2 += um[i] * z[i];
  }
  alfa = scal1 / scal2;
  for (int i = 0; i < Nuz; i++)
  {
    q[i] += alfa * z[i];
    r[i] -= alfa * um[i];
  }
  scal3 = 0;
  for (int i = 0; i < Nuz; i++)
  {
    scal3 += r[i] * r[i];
  }
  beta = scal3 / scal1;
  for (int i = 0; i < Nuz; i++)
  {
    z[i] = r[i] + beta * z[i];
  }
  residual = euclideanNorm(r) / bnorm;

  for (int k = 1; k < maxIter && residual > eps; k++)
  {
    scal1 = scal3;
    scal2 = 0;
    scal3 = 0;
    alfa = 0;
    beta = 0;
    mult(z, um);
    for (int i = 0; i < Nuz; i++)
    {
      scal2 += um[i] * z[i];
    }
    alfa = scal1 / scal2;
    for (int i = 0; i < Nuz; i++)
    {
      q[i] += alfa * z[i];
      r[i] -= alfa * um[i];
    }
    scal3 = 0;
    for (int i = 0; i < Nuz; i++)
    {
      scal3 += r[i] * r[i];
    }
    beta = scal3 / scal1;
    for (int i = 0; i < Nuz; i++)
    {
      z[i] = r[i] + beta * z[i];
    }
    residual = euclideanNorm(r) / bnorm;
  }
}

// Main function
int main()
{
  vector<nd> mesh;
  vector<el> elList;
  vector<double> time;
  double gamma;
  int fnum;
  input(mesh, elList, gamma, fnum, time);

  int Nuz = mesh.size();
  vector<vector<double>> M(3, vector<double>(3)), G(3, vector<double>(3));
  vector<double> locb(3, 0), F(Nuz, 0), q(Nuz, 0), q1(Nuz, 0), q2(Nuz, 0);
  vector<double> di(Nuz, 0), gg;
  vector<int> ig(Nuz + 1), jg;
  vector<vector<int>> list(Nuz);

  // Build sparse matrix portrait
  buildPortrait(list, ig, jg, elList);
  gg.resize(ig.back(), 0);

  // Initialize q2 and q1 (e.g., using exact solution or zero)
  for (int i = 0; i < Nuz; i++)
  {
    q2[i] = u(fnum, mesh[i]);
    q1[i] = u(fnum, mesh[i]);
  }

  // Time iteration
  for (int s = 2; s < time.size(); s++)
  {
    double t = time[s];
    double dt = time[s] - time[s - 1];

    // Reset global matrix and vector
    fill(di.begin(), di.end(), 0);
    fill(gg.begin(), gg.end(), 0);
    fill(F.begin(), F.end(), 0);

    // Assemble global system
    for (const auto &e : elList)
    {
      getM(M, gamma, e);
      getG(G, e);
      getb(locb, e, gamma, fnum, t, q1, q2, dt);

      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j <= i; j++)
        {
          double val = G[i][j] + (3.0 / (2.0 * dt)) * M[i][j];
          addElementToGlobal(e.nds[i].gl_num, e.nds[j].gl_num, val, di, gg, ig,
                             jg);
        }
        F[e.nds[i].gl_num] += locb[i];
      }
    }

    // Apply boundary conditions
    applyBC1(di, gg, F, ig, jg, mesh);

    // Solve system
    conjugateGradient(di, gg, ig, jg, F, q, 1e-10, 10000);

    // Update solutions
    q2 = q1;
    q1 = q;
    output(q, mesh, t);
  }

  return 0;
}