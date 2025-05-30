#include "header.h"

/*========================= I/O  FUNCTIONS =========================*/

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

/*==================== LOCAL MATRICES FUNCTIONS ====================*/

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
void getb(vector<double> &b, const el &e, double gamma, int s, double t,
          const vector<double> &q1, const vector<double> &q2, double dt)
{
  double det = detD(e);
  vector<double> p_loc(3);
  for (int i = 0; i < 3; i++)
  {
    p_loc[i] = f(s, e.nds[i], e); // Assuming f is defined as in your code
  }

  vector<vector<double>> M(3, vector<double>(3));
  getM(M, gamma, e);

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

/*=================== SPARSE MATRICES FUNCTIONS ===================*/
