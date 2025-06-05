#include "header.h"

void input(vector<nd> &mesh, vector<el> &elList, double &sigma, int &fnum)
{
  int n, k, node;
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
  int lambda;
  in >> lambda;
  for (int i = 0; i < n; i++)
    elList[i].lambda = lambda;
  in.close();

  // Read sigma
  in.open("../data/sigma.txt");
  in >> sigma;
  in.close();

  // Read fnum
  in.open("../data/fnum.txt");
  in >> fnum;
  in.close();
}

void input_timemesh(TimeMesh &time)
{
  int countTimeLayer = 0;

  ifstream timeMesh("../data/timeMesh.txt");

  timeMesh >> countTimeLayer;
  time.t.resize(countTimeLayer);

  for (int ti = 0; ti < countTimeLayer; ti++)
    timeMesh >> time.t[ti];

  timeMesh.close();
}

void input_split_timemesh(TimeMesh &time)
{
  auto &tNew = time.tNew;

  const int tSize = time.t.size();

  ifstream splittingTimeMesh("../data/timeMeshSplit.txt");

  int nk = 0;
  tNew.resize(1, time.t[0]);

  for (int ti = 0, j = 1; ti < tSize - 1; ti++, j++)
  {
    int countIntervals = 0;

    double coef = 0;
    double step = 0;

    splittingTimeMesh >> countIntervals >> coef;
    nk += countIntervals;
    tNew.resize(nk + 1);

    if (coef != 1)
    {
      double sumProgression = (pow(coef, countIntervals) - 1.) / (coef - 1.);
      step = (time.t[ti + 1] - time.t[ti]) / sumProgression;

      int jk = 1;
      for (j; j < nk; j++, jk++)
        tNew[j] = time.t[ti] + step * (pow(coef, jk) - 1.) / (coef - 1.);
    }
    else
    {
      step = (time.t[ti + 1] - time.t[ti]) / countIntervals;

      int jk = 1;
      for (j; j < nk; j++, jk++)
        tNew[j] = time.t[ti] + step * jk;
    }

    tNew[j] = time.t[ti + 1];
  }

  time.q.resize(time.tNew.size());
}

void input_boundary(vector<BoundaryConditions> &cond)
{
  int numEdgeConditions = 0;

  ifstream conditions("../data/boundaryConditions.txt");
  conditions >> numEdgeConditions;
  cond.resize(numEdgeConditions);

  for (int i = 0; i < numEdgeConditions; i++)
  {
    auto &VerticesNumbers = cond[i].VerticesNumbers;
    int numVertex = 0, type = 0, function = 0;

    conditions >> type >> function >> numVertex;

    cond[i].type = type;

    cond[i].function = function;
    VerticesNumbers.resize(numVertex);

    for (int k = 0; k < numVertex; k++)
      conditions >> VerticesNumbers[k];
  }
}

void print_solution(const TimeMesh &timemesh, const vector<nd> &mesh)
{
  int timeLayers = timemesh.q.size();
  int numNodes = mesh.size();

  for (int ti = 0; ti < timeLayers; ++ti)
  {
    cout << "=== Time step " << ti << ", t = " << timemesh.tNew[ti] << " ===\n";

    for (int i = 0; i < numNodes; ++i)
    {
      cout << "Node " << i << " (r = " << mesh[i].r << ", z = " << mesh[i].z
           << "): " << timemesh.q[ti][i] << '\n';
    }

    cout << '\n';
  }
}

void out_solution(TimeMesh &timemesh, vector<nd> &mesh, int flag)
{
  ofstream out("../results/3t/conv_t^4.txt");
  // ofstream out("../results/3t/conv_cos.txt");
  out << std::scientific << std::setprecision(8);

  int timeLayers = timemesh.q.size();
  int numNodes = mesh.size();

  for (int ti = 0; ti < timeLayers; ++ti)
  {
    out << "t = " << timemesh.tNew[ti] << "\n";
    out << "       r                 z               t               u     "
           " "
           "      "
           "   u*              |u-u*|\n";

    for (int i = 0; i < numNodes; ++i)
    {
      double u_num = timemesh.q[ti][i]; // Численное решение
      double u_analytic =
          u(mesh[i], timemesh.tNew[ti], flag);     // Аналитическое решение
      double error = std::abs(u_num - u_analytic); // Погрешность

      out << std::setw(16) << mesh[i].r << std::setw(17) << mesh[i].z
          << std::setw(16) << timemesh.tNew[ti] << std::setw(17) << u_num
          << std::setw(16) << u_analytic << std::setw(17) << error << "\n";
    }
    out << "\n";
  }

  out.close();
}