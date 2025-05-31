# Math Physics Equations Course Project   
## Task   
### FEM for parabolic two-dimensional equation with three-layer implicit time scheme in cylindric coordinates on triangles. Linear basis functions.   

## To-do list
- [x] Check initial code
- [x] make input function
- [ ] Add sparse matrix format
- [ ] Add building of global matrix & vector in sparse format
- [ ] Add sparse matrix solver (CGM)
- [ ]  Add implicit 3-layer time scheme
- [ ]  Add 1st boundary conditionals
- [ ]  Add 2nd boundary conditionals
- [ ]  Rewrite ```main()``` so program will work
- [ ]  Refactor the code
- [ ]  Run the tests (maybe automate testing)
- [ ]  Make report

## Input files
`data/els.txt` - elements   
`data/nds.txt` - nodes of the mesh   
`data/fnum.txt` - number of the test function _f_   
`data/lambda.txt` - _lambda_ on each finite element   
`data/gamma.txt` - _gamma_    
`data/c1.txt` - first boundary conditions   
`data/timeMesh.txt` - **Временные слои (время). В первой строке целое число k – количество временных слоев, далее k значений временных слоев.**   
`data/timeCoef.txt` - **Разбиение временной сетки. В нем идут k-1 пар чисел (для каждого слоя), первое из которых значит на сколько интервалов нужно разбить интервал времени, второе – коэффициент разбиения (растяжения сжатия).**   
`data/q0.txt` - _q0_ initial vector???   

## Structures and data in program
### Structures
**`nd`:**   
`int gl_num` - global number of the node   
`double r, z` - coordinates of the node (r,z)   

**`el`:**   
`double lambda` - lambda value on the element   
`nd nds[3]` - nodes of the element   

**`timelayer`:**   
`double t` -  value of current timestep   
`int intervals` - number of intervals on current timestep   
`double coef` - coefficient of intervals on current time step   
   
### Data
`vector<nd> mesh` - _n_ global nodes   
`vector<el> elList` - elements   
`vector<timelayer> time` - _k_ time layers   
`double gamma` - _gamma_ (const)    
`int fnum` - number of the test function _f_    


## Files
`main.cc` - main file of the program
## Functions
to be continued...

## Algorithm
1. input
2. portrait
3. calculating q0 and q1
4. going through all other time layers
5. on each time layer new global matrix is being made
6.  solve SLAE
7.  save result to the array

После считывания входных файлов происходит построение портрета глобальной матрицы по информации из сетки (примерно описано выше). Далее запускается функция, решающая задачу, в которой сначала вычисляются значение весов на первых 2-х слоях, далее начинаем идти по временным слоям. Для каждого временного слоя собирается глобальная матрица и вектор правой части с помощью добавления локальных матриц и вектора правой части на каждом элементе в глобальные матрицу и вектор правой части и учета краевых условий. После сборки глобальной матрицы и вектора правой части, решается СЛАУ с помощью метода ЛОС с LU предобуславливанием. Полученный вектор решения СЛАУ сохраняется в массив, содержащий массивы весов разложения с каждого временного слоя.