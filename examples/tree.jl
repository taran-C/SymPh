using SymbolicPhysics.Arrays

x = ArrayVariable("x")
a = 2*x + x*x
b = 3*x
y = ArrayVariable("a") * ArrayVariable("b") #add names to Expression and this will go away as y = a+b AS IT SHOULD

#----To be done automatically with to_tree--- 
#Different nodes to represent x mean different places that depend on it, but shouldn't be a problem (a couple of nodes most importantly represents a dependecy)
xa = DepNode("x", x)
xb = DepNode("x", x)

an = DepNode("a", a)
bn = DepNode("b", b)

yn = DepNode("y", y)

addchild!(an, xa)
addchild!(bn, xb)

addchild!(yn, an)
addchild!(yn, bn)
#--------------------------------------------

seq = to_sequence!(yn)
println(string(seq))

func!, funcstr = Arrays.to_kernel(seq)
println(funcstr)

#Mesh
nx = 10
ny = 10
nh = 2

#Input Variables
x = ones(nx,ny)

#Diags
b = zeros(nx, ny)
a = zeros(nx, ny)
y = zeros(nx, ny)

display(y)

func!(;nx=nx, ny=ny, nh=nh, x=x, a=a, b=b, y=y)

display(y)
