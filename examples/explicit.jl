using SymbolicPhysics.Maths
import SymbolicPhysics.Arrays

#Defining our equation
U = VectorVariable{Dual}("U")
a = FormVariable{2,Primal}("a")
dta = ExteriorDerivative("dta", InteriorProduct(U,a))
print(dta)
#-----This part will be encapsulated into a function that automatically optimizes ----------
#Transforming the Forms expression into an Expression on arrays
exprs = explicit(dta)
println("Developped expression :")
println(string(exprs))

#Transforming our Expression into a dependency tree
tree = Arrays.to_deptree!(Set{String}(["ι_U_a_x","ι_U_a_y"]), exprs)
println("Tree view")
println(string(tree))
println("Graphviz view of tree")
println(Arrays.to_graphviz(tree))

#Transforming our dependency tree into a sequence of expressions to compute
seq = Arrays.to_sequence!(tree)
println("Corresponding Sequence :")
println(string(seq)*"\n")

#Generating the final function
func!, funcstr = Arrays.to_kernel(seq)
println("Generated code :")
println(funcstr)
#-------------------------------------------------------------------------------------------

#Testing the function
nx = 14
ny = 14
nh = 2

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

mesh = Arrays.Mesh(nx, ny, nh, msk, 10, 10)

A = mesh.xc .+ mesh.yc
A .*= msk

U_X = deepcopy(mesh.xc)
U_X .*= mesh.mskx
U_Y = deepcopy(mesh.yc)
U_Y .*= mesh.msky

DTA = zeros(nx, ny)

ι_U_a_x = zeros(nx, ny)
ι_U_a_y = zeros(nx, ny)

func!(;nx=nx, ny=ny, nh=nh, a=A, U_X=U_X, U_Y=U_Y, dta=DTA, ι_U_a_x = ι_U_a_x, ι_U_a_y = ι_U_a_y, mskx = mesh.mskx, msky = mesh.msky, mskv = mesh.mskv)

display(A)
display(DTA)
