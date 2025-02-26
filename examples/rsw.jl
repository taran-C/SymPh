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
@Let f = FormVariable{2, Dual}() #Coriolis (f* so times A)

#Time derivative
@Let du = -InteriorProduct(U, zeta + f) - ExteriorDerivative(p + k) #du = -i(U, ζ* + f*) - d(p + k)
@Let dh = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))

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
tree = Arrays.to_deptree!(Set{String}(["zeta"]), exprs)
#println("Tree view")
#println(string(tree))
#println("Graphviz view of tree")
#println(Arrays.to_graphviz(tree))

#Transforming our dependency tree into a sequence of expressions to compute
seq = Arrays.to_sequence!(tree)
#println("Corresponding Sequence :")
#println(string(seq)*"\n")

#Generating the final function
rsw!, rswstr = Arrays.to_kernel(seq)
#println("Generated code :")
println(rswstr)
#-------------------------------------------------------------------------------------------

#Testing the function

#Defining the Mesh
nx = 100
ny = 100
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

Lx, Ly = (1,1)
mesh = Arrays.Mesh(nx, ny, nh, msk, Lx, Ly)

#Initial Conditions
h0 = 0.15
H = 1
sigma = 0.05
d = 1.4*sigma 
gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
h = zeros((nx,ny))

config = "vortex"

for i in nh+1:nx-nh, j in nh+1:ny-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]

	if config == "dipole"
		h[i,j] = (H + h0 * (gaussian(x, y, 0.45, 0.5, sigma) - gaussian(x,y, 0.55, 0.5, sigma))) * mesh.A[5,5]
	elseif config == "vortex"
		h[i,j] = (H + h0 * gaussian(x, y, 0.5, 0.5, sigma)) * mesh.A[5,5]
	end
end

f = 5 .* ones((nx,ny)) .* mesh.A .* mesh.msk2d

u_x = zeros((nx, ny))
u_y = zeros((nx, ny))

zeta = zeros((nx,ny))

dh1 = zeros((nx, ny))
dh2 = zeros((nx, ny))
dh3 = zeros((nx, ny))
du_x1 = zeros((nx, ny))
du_x2 = zeros((nx, ny))
du_x3 = zeros((nx, ny))
du_y1 = zeros((nx, ny))
du_y2 = zeros((nx, ny))
du_y3 = zeros((nx, ny))

#Integrator TODO put somewhere else or use an externalized one
function rk3step!(dt, mesh, f, h, u_x, u_y, zeta,
		dh1, dh2, dh3,
		du_x1, du_x2, du_x3,
		du_y1, du_y2, du_y3)
	
	rsw!(mesh ;f=f, h=h, u_x=u_x, u_y=u_y, zeta=zeta, dh=dh1, du_x=du_x1, du_y=du_y1)
	h .+= dt .* dh1
	u_x .+= dt .* du_x1
	u_y .+= dt .* du_y1

	rsw!(mesh ;f=f, h=h, u_x=u_x, u_y=u_y, zeta=zeta, dh=dh2, du_x=du_x2, du_y=du_y2)
	h .+= dt .* (-3/4 .* dh1 + 1/4 .* dh2)
	u_x .+= dt .* (-3/4 .* du_x1 + 1/4 .* du_x2)
	u_y .+= dt .* (-3/4 .* du_y1 + 1/4 .* du_y2)

	rsw!(mesh ;f=f, h=h, u_x=u_x, u_y=u_y, zeta=zeta, dh=dh3, du_x=du_x3, du_y=du_y3)
	h .+= dt .* (-1/12 .* dh1 - 1/12 .* dh2 + 2/3 .* dh3)
	u_x .+= dt .* (-1/12 .* du_x1 - 1/12 .* du_x2 + 2/3 .* du_x3)
	u_y .+= dt .* (-1/12 .* du_y1 - 1/12 .* du_y2 + 2/3 .* du_y3)
end

#Plotting
plot = true
var_to_plot = h

if plot
	using GLMakie
	global fig = Figure(size = (600, 600))
	global ax = Axis(fig[1,1])
	global q = Observable(var_to_plot[nh+1:nx-nh, nh+1:ny-nh])
	global hm = heatmap!(ax, q)#, colorrange = (0.0001, 0.0002))
	display(fig)
end

#TimeLoop
tend = 10
maxite = 3000
ite = 0
dt = 0.15 * minimum(mesh.dx[nh+1:nx-nh, nh+1:ny-nh])

tstart = time()
for t in 0:dt:tend
	if ite >= maxite
		break
	end

	#Actual progress
	rk3step!(dt, mesh, f, h, u_x, u_y, zeta,
		 dh1, dh2, dh3,
		 du_x1, du_x2, du_x3,
		 du_y1, du_y2, du_y3)

	if ite % 10 == 0
		println("ite : $(ite)/$(maxite), t : $(t)/$(tend)")
		if plot
			q[] = var_to_plot[nh+1:nx-nh, nh+1:ny-nh]
		end
		sleep(0.1)
	end
	if ite%100 == 0
		save("rsw_$(ite).png", fig)
	end

	global ite += 1
end
println(time()-tstart)
