# Math Physics Equations Course Project   
## Task   
### FEM for parabolic two-dimensional equation with three-layer implicit time scheme in cylindric coordinates on triangles. Linear basis functions.   

## To-do list
- [x] Check initial code
- [x] make input function
- [x] Add sparse matrix format
- [x] Add building of global matrix & vector in sparse format
- [x] Add sparse matrix solver (LOS)
- [x] Add implicit 3-layer time scheme
- [x] Add 1st boundary conditionals
- [ ] Add 2nd boundary conditionals
- [x] Rewrite ```main()``` so program will work
- [ ] Test on different functions
- [ ] Find out why there is such big error in the inner mode
- [ ] Refactor the code
- [ ] Run the tests (maybe automate testing)
- [ ] Save the results of testing
- [ ] Make report

## Input files
`data/els.txt` - elements   
`data/nds.txt` - nodes of the mesh   
`data/fnum.txt` - number of the test function _f_   
`data/lambda.txt` - _lambda_ on each finite element   
`data/gamma.txt` - _gamma_    
`data/boundaryConditions.txt` - first boundary conditions   
`data/timeMesh.txt` -  
`data/timeCoef.txt` - 

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
