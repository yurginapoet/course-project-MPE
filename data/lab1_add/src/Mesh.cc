#include "Mesh.h"

void Mesh::Clear()
{
  nodes.clear();
  first.clear();
  second.clear();
  b.clear();
  u.clear();
  fictnum.clear();
  matrix.clear();

  num = 0;
  nx = 0;
  ny = 0;
  m = 0;
  hx = 0.0;
  hy = 0.0;
  kx = 0.0;
  ky = 0.0;
  rx = 0.0;
  ry = 0.0;
  lambda = 0.0;
  gamma = 0.0;
  eps = 0.0;
  w = 0.0;
  iternum = 0;
  testnum = 0;

  area.xmin = 0.0;
  area.xmax = 0.0;
  area.ymin = 0.0;
  area.ymax = 0.0;
}

// input the data
void Mesh::Input(string fin, string flg, string fparam)
{
  ifstream in("data/blocks.txt");
  in >> N;
  blocks.resize(N);
  blocks[0].resize(4);
  in >> blocks[0][0] >> blocks[0][1] >> blocks[0][2] >> blocks[0][3];

  for (int i = 1; i < N; i++)
  {
    blocks[i].resize(4);
    in >> blocks[i][0] >> blocks[i][1] >> blocks[i][2] >> blocks[i][3];
  }
  in.close();

  in.open(fin);

  in >> nx >> ny;
  in >> kx >> ky;
  in >> area.xmin >> area.xmax >> area.ymin >> area.ymax;
  // in >> rx >> ry;
  in.close();
  in.open(flg);
  in >> lambda >> gamma;
  in.close();
  in.open(fparam);
  in >> eps >> w >> iternum;
  in.close();
}

// calculate step by x and y
void Mesh::CalcStep()
{
  if (kx == 1)
    hx = (area.xmax - area.xmin) / nx;
  else
    hx = Step(nx, kx, area.xmax, area.xmin);

  if (ky == 1)
    hy = (area.ymax - area.ymin) / ny;
  else
    hy = Step(ny, ky, area.ymax, area.ymin);
}

// make mesh of nodes
void Mesh::Nodes()
{
  num = (nx + 1) * (ny + 1);
  nodes.resize(num);

  for (int k = 0; k < num; k++)
  {
    int row = k / (nx + 1);
    int col = k - row * (nx + 1);

    d xnode = 0, ynode = 0;

    for (int i = 0; i < col; i++)
      xnode += hx * pow(kx, i);
    for (int i = 0; i < row; i++)
      ynode += hy * pow(ky, i);

    if (col == nx)
      nodes[k].x = area.xmax; // right
    else
      nodes[k].x = area.xmin + xnode; // inner by x

    if (row == ny)
      nodes[k].y = area.ymax; // upper
    else
      nodes[k].y = area.ymin + ynode; // inner by y
  }
}

// EDIT!!!!
void Mesh::FindFictive()
{
  for (int i = 0; i < nodes.size(); i++)
  {
    int f = 0;
    for (int j = 0; j < N; j++)
      if (nodes[i].x <= blocks[j][1] && nodes[i].x >= blocks[j][0] &&
          nodes[i].y <= blocks[j][3] && nodes[i].y >= blocks[j][2])
      {
        f = 1;
        break;
      }
    if (!f)
    {
      fictnum.resize(fictnum.size() + 1);
      fictnum.back() = i;
    }
  }
}

// check if node is fictive
bool Mesh::IsFictive(int k)
{
  if (binary_search(fictnum.begin(), fictnum.end(), k))
    return true;
  else
    return false;
}

// check if node is inner
bool Mesh::IsInnerNode(int i, int j, int k)
{
  if (i != 0 && i != ny && j != 0 && j != nx && !IsFictive(k + 1) &&
      !IsFictive(k - 1) && !IsFictive(k + m) && !IsFictive(k - m) &&
      !IsFictive(k + m + 1) && !IsFictive(k + m - 1) && !IsFictive(k - m + 1) &&
      !IsFictive(k - m - 1))
    return true;
  else
    return false;
}

// check if node belongs to the first conditional
bool Mesh::IsFirst(int k)
{
  if (binary_search(first.begin(), first.end(), k))
    return true;
  else
    return false;
}

// check if node belongs to the second conditional
bool Mesh::IsSecond(int k)
{
  if (binary_search(second.begin(), second.end(), k))
    return true;
  else
    return false;
}

// input and add 1st and 2nd boundary conditions
void Mesh::Boundary(string ffirst, string fsecond)
{
  int number; // number of blocks with conditional
  int side;   // border

  ifstream infile(fsecond);
  infile >> number;

  for (int i = 0; i < number; i++)
  {
    double f, s; // Границы отрезка с заданным к. у.
    int bl, gr;  // Блок(нумерация с 1) и грань на которой оно задано
    infile >> bl >> gr;
    switch (gr)
    {
      // down
    case 1:
      for (int j = 0; j < nodes.size(); j++)
        if (nodes[j].x >= blocks[bl - 1][0] &&
            nodes[j].x <= blocks[bl - 1][1] && nodes[j].y == blocks[bl - 1][2])
        {
          second.resize(second.size() + 1);
          second.back() = j;
        }
      break;
      // right
    case 2:
      for (int j = 0; j < nodes.size(); j++)
        if (nodes[j].y >= blocks[bl - 1][2] &&
            nodes[j].y <= blocks[bl - 1][3] && nodes[j].x == blocks[bl - 1][1])
        {
          second.resize(second.size() + 1);
          second.back() = j;
        }
      break;
      // up
    case 3:
      for (int j = 0; j < nodes.size(); j++)
        if (nodes[j].x >= blocks[bl - 1][0] &&
            nodes[j].x <= blocks[bl - 1][1] && nodes[j].y == blocks[bl - 1][3])
        {
          second.resize(second.size() + 1);
          second.back() = j;
        }
      break;
      // left
    case 4:
      for (int j = 0; j < nodes.size(); j++)
        if (nodes[j].y >= blocks[bl - 1][2] &&
            nodes[j].y <= blocks[bl - 1][3] && nodes[j].x == blocks[bl - 1][0])
        {
          second.resize(second.size() + 1);
          second.back() = j;
        }
      break;
    }
  }

  infile.close();
  set<d> temp(second.begin(), second.end());
  second.assign(temp.begin(), temp.end());
  temp.clear();

  infile.open(ffirst);
  infile >> number;

  for (int i = 0; i < number; i++)
  {
    double f, s; // Границы отрезка с заданным к. у.
    int bl, gr;  // Блок(нумерация с 1) и грань на которой оно задано
    infile >> bl >> gr;
    switch (gr)
    {
    case 1:
      for (int j = 0; j < nodes.size(); j++)
        if (nodes[j].x >= blocks[bl - 1][0] &&
            nodes[j].x <= blocks[bl - 1][1] && nodes[j].y == blocks[bl - 1][2])
        {
          first.resize(first.size() + 1);
          first.back() = j;
        }
      break;
    case 2:
      for (int j = 0; j < nodes.size(); j++)
        if (nodes[j].y >= blocks[bl - 1][2] &&
            nodes[j].y <= blocks[bl - 1][3] && nodes[j].x == blocks[bl - 1][1])
        {
          first.resize(first.size() + 1);
          first.back() = j;
        }
      break;
    case 3:
      for (int j = 0; j < nodes.size(); j++)
        if (nodes[j].x >= blocks[bl - 1][0] &&
            nodes[j].x <= blocks[bl - 1][1] && nodes[j].y == blocks[bl - 1][3])
        {
          first.resize(first.size() + 1);
          first.back() = j;
        }
      break;
    case 4:
      for (int j = 0; j < nodes.size(); j++)
        if (nodes[j].y >= blocks[bl - 1][2] &&
            nodes[j].y <= blocks[bl - 1][3] && nodes[j].x == blocks[bl - 1][0])
        {
          first.resize(first.size() + 1);
          first.back() = j;
        }
      break;
    }
  }

  infile.close();

  temp.insert(first.begin(), first.end());
  first.assign(temp.begin(), temp.end());
  temp.clear();
}
// remove 2nd conditional which belong to the 1st conditional
void Mesh::RemoveSecondConditionals()
{
  int i = 0, j = 0;
  while (i < second.size() && j < first.size())
  {
    if (second[i] < first[j])
      i++;
    else if (second[i] > first[j])
      j++;
    else
      second.erase(second.begin() + i);
  }
}

// fill in the matrix
void Mesh::FillMatrix()
{
  m = nx + 1;
  matrix.resize(5);

  matrix[2].resize(num, 0);
  matrix[1].resize(num - 1, 0);
  matrix[3].resize(num - 1, 0);
  matrix[0].resize(num - m, 0);
  matrix[4].resize(num - m, 0);

  b.resize(num, 0);
  u.resize(num, 0);

  vector<d> Hx(nx), Hy(ny);

  Hx[0] = hx;
  Hy[0] = hy;
  for (int i = 1; i < nx; i++)
    Hx[i] = Hx[i - 1] * kx;
  for (int i = 1; i < ny; i++)
    Hy[i] = Hy[i - 1] * ky;

  for (int i = 0; i <= ny; i++)
  {
    for (int j = 0; j <= nx; j++)
    {
      int k = i * m + j;
      if (IsFictive(k))
      {
        matrix[2][k] = 1;
        b[k] = 0;
      }
      else if (IsInnerNode(i, j, k))
      {
        matrix[2][k] =
            lambda * (2 / (Hx[j] * Hx[j - 1]) + 2 / (Hy[i] * Hy[i - 1])) +
            gamma;
        // if (j > 0)
        matrix[1][k - 1] = -2 * lambda / (Hx[j - 1] * (Hx[j] + Hx[j - 1]));
        // if (j < nx)
        matrix[3][k] = -2 * lambda / (Hx[j] * (Hx[j] + Hx[j - 1]));
        // if (i > 0)
        matrix[0][k - m] = -2 * lambda / (Hy[i - 1] * (Hy[i] + Hy[i - 1]));
        // if (i < ny)
        matrix[4][k] = -2 * lambda / (Hy[i] * (Hy[i] + Hy[i - 1]));
        b[k] = F(nodes[k].x, nodes[k].y, testnum);
      }
      else
      {
        if (IsFirst(k))
        {
          matrix[2][k] = 1;
          b[k] = U(nodes[k].x, nodes[k].y, testnum);
        }
        else if (IsSecond(k))
        {
          if (i == 0 || IsFictive(k - m))
          {
            matrix[2][k] = -lambda / Hy[i];
            matrix[4][k] = lambda / Hy[i];
            b[k] = -Theta(nodes[k].x, nodes[k].y, testnum);
          }
          else if (i == ny || IsFictive(k + m))
          {
            matrix[2][k] = lambda / Hy[i - 1];
            matrix[4][k - m] = -lambda / Hy[i - 1];
            b[k] = Theta(nodes[k].x, nodes[k].y, testnum);
          }
          else if (j == 0 || IsFictive(k - 1))
          {
            matrix[2][k] = lambda / Hx[j];
            matrix[3][k] = -lambda / Hx[j];
            b[k] = Theta(nodes[k].x, nodes[k].y, testnum);
          }
          else if (j == nx || IsFictive(k + 1))
          {
            matrix[2][k] = lambda / Hx[j - 1];
            matrix[1][k - 1] = -lambda / Hx[j - 1];
            b[k] = Theta(nodes[k].x, nodes[k].y, testnum);
          }
        }
      }
    }
  }
}

// print all information from class
void Mesh::Out()
{
  cout << "Mesh variables:" << endl;

  // Вывод простых переменных
  cout << "nx: " << nx << endl;
  cout << "ny: " << ny << endl;
  cout << "m: " << m << endl;
  cout << "hx: " << hx << endl;
  cout << "hy: " << hy << endl;
  cout << "kx: " << kx << endl;
  cout << "ky: " << ky << endl;
  cout << "lambda: " << lambda << endl;
  cout << "gamma: " << gamma << endl;
  cout << "eps: " << eps << endl;
  cout << "w: " << w << endl;
  cout << "iternum: " << iternum << endl;

  // Вывод вектора nodes
  cout << "nodes:" << endl;
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    cout << "  Node " << i << ": (" << nodes[i].x << ", " << nodes[i].y << ")"
         << endl;
  }

  // Вывод вектора rects
  cout << "rx & ry:" << endl;
  cout << rx << ' ' << ry << endl;
  cout << defaultfloat;
  // Вывод векторов first, second, b, u
  auto PrintVector = [](const string &name, const vector<d> &vec)
  {
    cout << name << ":" << endl;
    for (size_t i = 0; i < vec.size(); ++i)
    {
      cout << "  " << name << "[" << i << "]: " << vec[i] << endl;
    }
  };

  auto PrintV = [](const string &name, const vector<int> &vec)
  {
    cout << name << ":" << endl;
    for (size_t i = 0; i < vec.size(); ++i)
    {
      cout << "  " << name << "[" << i << "]: " << vec[i] << endl;
    }
  };

  PrintV("first", first);
  PrintV("second", second);
  cout << scientific << setprecision(3);
  PrintVector("b", b);
  PrintVector("u", u);

  // Вывод вектора fictnum
  cout << "fictnum:" << endl;
  for (size_t i = 0; i < fictnum.size(); ++i)
  {
    cout << "  fictnum[" << i << "]: " << fictnum[i] << endl;
  }

  // Вывод area
  cout << "area: xmin=" << area.xmin << ", xmax=" << area.xmax
       << ", ymin=" << area.ymin << ", ymax=" << area.ymax << endl;

  // Вывод матрицы matrix
  cout << "matrix:" << endl;
  cout << defaultfloat;
  for (size_t i = 0; i < matrix.size(); ++i)
  {
    cout << "  Row " << i << ": ";
    for (size_t j = 0; j < matrix[i].size(); ++j)
    {
      cout << matrix[i][j] << " ";
    }
    cout << endl;
  }
}

// Hauss-Zeidel solver
void Mesh::HaussZeidel(string filename)
{
  ofstream o(filename);
  vector<d> Ax;
  int num = nodes.size();
  int k = 0;
  Ax.resize(num);
  do
  {
    for (int i = 0; i < num; i++)
      Ax[i] = b[i] - Iteration(i);
    k++;
    o << "iteration: " << k << '\t' << "residual: " << scientific
      << setprecision(15) << Norm(num, Ax) / Norm(num, b) << endl;
  } while (k < iternum && Norm(num, Ax) / Norm(num, b) >= eps);
}

// iterational process
d Mesh::Iteration(int i)
{
  d sum = 0, a = 0;
  for (int j = 0; j < num; j++)
  {
    int r = j - i;
    if (r == 0)
      a = matrix[2][i];
    else if (r == 1)
      a = matrix[3][i];
    else if (r == m)
      a = matrix[4][i];
    else if (r == -1)
      a = matrix[1][j];
    else if (r == -m)
      a = matrix[0][j];
    else
      a = 0;
    sum += a * u[j];
  }
  u[i] = u[i] + w / matrix[2][i] * (b[i] - sum);
  return sum;
}

// print results of executing the program
void Mesh::Results(string filename)
{
  ofstream o(filename);
  if (!o)
  {
    cerr << "Error opening file: " << filename << endl;
    return;
  }

  o << setw(14) << "x" << setw(14) << "y" << setw(26) << "u" << setw(26) << "u*"
    << setw(26) << "|u* - u|" << endl;

  for (int i = 0; i < nodes.size(); i++)
  {
    o << scientific << setprecision(4) << setw(14) << nodes[i].x << setw(14)
      << nodes[i].y;

    if (binary_search(fictnum.begin(), fictnum.end(), i))
      o << setw(26) << "-" << setw(26) << "-" << setw(26) << "-";
    else
    {
      double u_star = U(nodes[i].x, nodes[i].y, testnum);
      double discrepancy = abs(u_star - u[i]);

      o << scientific << setprecision(15) << setw(26) << u[i] << setw(26)
        << u_star << setw(26) << discrepancy;
    }

    o << endl;
  }

  o.close();
}