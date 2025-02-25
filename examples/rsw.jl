import SymbolicPhysics: @Let
using SymbolicPhysics.Maths
import SymbolicPhysics.Arrays

#Defining our equation
@Let h = FormVariable{2, Primal}() #Height * A (h* technically)
@Let u = FormVariable{1, Dual}() #Transported velocity

@Let U = Sharp(u) # U = u#
@Let k = 0.5 * Hodge(InnerProduct(u,u)) #k = 0.5 * hodge(innerproduct(u,u))
@Let p = Hodge(h) # p = *(g(h*+b*))
@Let zeta = ExteriorDerivative(u) # ζ* = du
@Let f = FormVariable{2, Dual}() #Coriolis

#Time der
@Let du = -InteriorProduct(U, zeta + f) - ExteriorDerivative(p + k) #du = -i(U, ζ* + f*) - d(p + k)
@Let dh = ExteriorDerivative(InteriorProduct(U, h)) #dh = Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))

#Checking the typings
@assert InnerProduct(u,u) isa Form{2, Primal}
@assert U isa Vect{Dual}
@assert k isa Form{0, Dual}
@assert p isa Form{0, Dual}

#-----This part will be encapsulated into a function that automatically optimizes ----------
#Transforming the Forms expression into an Expression on arrays
exprs = [explicit(du); explicit(dh)]
#println("Developped expression :")
#println(string(exprs))

#Transforming our Expression into a dependency tree
tree = Arrays.to_deptree!(Set{String}([]), exprs)
#println("Tree view")
#println(string(tree))
#println("Graphviz view of tree")
#println(Arrays.to_graphviz(tree))

#Transforming our dependency tree into a sequence of expressions to compute
seq = Arrays.to_sequence!(tree)
#println("Corresponding Sequence :")
#println(string(seq)*"\n")

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
plot = false

if plot :
	using GLMakie
	fig = Figure(size = (600, 600))
	ax = Axis(fig[1,1])
	q = Observable(A)

#TimeLoop
tend = 50
dt = 0.01

tstart = time()
for t in 0:dt:tend
	func!(mesh ;a=A, U_X=U_X, U_Y=U_Y, dta=DTA, ι_U_a_x = ι_U_a_x, ι_U_a_y = ι_U_a_y)
	#func!(mesh ;a=A, U_X=U_X, U_Y=U_Y, dta=DTA)
	A .-= dt .* DTA
end
println(time()-tstart)

#display(A)
#display(DTA)
