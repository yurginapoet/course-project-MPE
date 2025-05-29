#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

// Структура для представления точки
struct nd
{
  double x, y;
};

// Структура для представления треугольника
struct el
{
  int v1, v2, v3; // Индексы вершин треугольника
};

// Функция для генерации регулярной сетки треугольников
void generate_mesh(int nx, int ny, vector<nd> &nds, vector<el> &els)
{
  // Генерация узлов сетки
  for (int j = 0; j <= ny; ++j)
  {
    for (int i = 0; i <= nx; ++i)
    {
      nds.push_back({1.0 + static_cast<double>(i) * 4.0 / nx,
                     1.0 + static_cast<double>(j) * 4.0 / ny});
    }
  }

  // Генерация треугольников
  for (int j = 0; j < ny; ++j)
  {
    for (int i = 0; i < nx; ++i)
    {
      int base = j * (nx + 1) + i;
      els.push_back({base, base + 1, base + nx + 1});
      els.push_back({base + 1, base + nx + 2, base + nx + 1});
    }
  }
}

// Функция для проверки, является ли узел граничным
bool is_boundary_node(nd &nd, int nx, int ny)
{
  return nd.x == 1.0 || nd.x == 5.0 || nd.y == 1.0 || nd.y == 5.0;
}

int main()
{
  int nx, ny;
  int c1, g, l, fnum;
  vector<nd> nds;
  vector<el> els;
  vector<int> bnds;
  cout << "nx & ny: " << endl;
  cin >> nx >> ny;
  cout << "gamma, lambda, fnum, c1: " << endl;
  cin >> g >> l >> fnum >> c1;

  // Генерация сетки
  generate_mesh(nx, ny, nds, els);

  ofstream f("../data/nds.txt");

  // Вывод узлов
  f << nds.size() << endl;
  for (const auto &nd : nds)
    f << nd.x << ' ' << nd.y << endl;

  f.close();

  f.open("../data/c1.txt");
  // c1
  for (int i = 0; i < nds.size(); ++i)
    if (is_boundary_node(nds[i], nx, ny))
      bnds.push_back(i);

  int n = bnds.size() - 1;

  f << n + 1 << endl;
  bool flag = false;
  for (int i = 0; i < n; i++)
  {
    if (i < nx || flag)
      f << bnds[i] << ' ' << bnds[i + 1] << ' ' << c1 << endl;
    else if (i == n - nx - 1)
    {
      f << bnds[i] << ' ' << bnds[n] << ' ' << c1 << endl;
      flag = true;
    }
    else
      f << bnds[i] << ' ' << bnds[i + 2] << ' ' << c1 << endl;
  }
  f << bnds[0] << ' ' << bnds[nx + 1] << ' ' << c1 << endl;

  f.close();

  f.open("../data/els.txt");
  // Вывод треугольников
  f << els.size() << endl;
  for (const auto &el : els)
    f << el.v1 << " " << el.v2 << " " << el.v3 << endl;
  f.close();

  f.open("../data/gamma.txt");
  f << g;

  f.close();
  f.open("../data/lambda.txt");
  for (int i = 0; i < els.size(); i++)
    f << l << ' ';

  f.close();
  f.open("../data/fnum.txt");
  f << fnum;
  f.close();

  return 0;
}