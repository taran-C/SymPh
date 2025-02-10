using SymbolicPhysics.Maths
import SymbolicPhysics.Arrays

#using GLMakie

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
#tree = Arrays.to_deptree!(Set{String}([]), exprs)
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
nx = 100
ny = 100
nh = 2

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

Lx, Ly = (10,10)
mesh = Arrays.Mesh(nx, ny, nh, msk, 10, 10)

#Grid
A = mesh.mskv .* (mod.(floor.(mesh.xc ./ (6 * mesh.dx-1)), 2) + mod.(floor.(mesh.yc ./ (6 * mesh.dy-1)), 2) .- 1)

U_X = cos.(pi * (mesh.xc ./Lx .- 0.5)) .* sin.(pi * (mesh.yc ./Ly .- 0.5)) .* mesh.mskx
U_Y = -cos.(pi * (mesh.yc ./Ly .- 0.5)) .* sin.(pi * (mesh.xc ./Lx .- 0.5)) .* mesh.msky

DTA = zeros(nx, ny)

ι_U_a_x = zeros(nx, ny)
ι_U_a_y = zeros(nx, ny)


#Heatmap
#fig = Figure(size = (600, 600))
#ax = Axis(fig[1,1])

#q = Observable(A)
#h = heatmap!(ax, q, colorrange=(-1.5,1.5))
#Colorbar(fig[1,2], h)
#display(fig)

#TimeLoop
tend = 50
dt = 0.01

tstart = time()
for t in 0:dt:tend
	func!(mesh ;a=A, U_X=U_X, U_Y=U_Y, dta=DTA, ι_U_a_x = ι_U_a_x, ι_U_a_y = ι_U_a_y)
	#func!(mesh ;a=A, U_X=U_X, U_Y=U_Y, dta=DTA)
	#A .-= dt .* DTA
	#q[] = A
	#sleep(0.00001)
end
println(time()-tstart)

#display(A)
#display(DTA)
