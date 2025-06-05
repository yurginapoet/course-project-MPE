# Method of finite differences (basic)   

Information:

## Input files:
### **`"data/in.txt"`**   
Information about area.   
`nx`, `ny` - number of intervals by x, y;  
`kx`, `ky` - coefficient of compression by x, y;  
`x_min x_max y_min y_max` - area boundary;   
`rx ry` - coordinates of point defining area.   

### **`"data/lg.txt"`**  
Information about parameters lambda and gamma.   
`lambda` `gamma` - parameters  

### **`"data/param.txt"`**   
Information about parameters for Hauss-Zeidel solver.   
`eps` - epsilon   
`w` - relaxation parameter   
`iternum` - max number of iterations   

### **`"data/first.txt"`**   
Information about first boundary conditions.   
`N` - number of boundaries with 1st conditional.      
N lines `side`, where `side` is 0-5;    

### **`"data/second.txt"`**   
Information about first boundary conditions.   
`N` - number of boundaries with 2nd conditional.      
N lines `side`, where `side` is 0-5;    