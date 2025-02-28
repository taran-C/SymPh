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
@Let pv = (f + zeta) / h #TODO check what pv should be

#Time derivative
@Let du = -InteriorProduct(U, zeta + f) - ExteriorDerivative(p + k) #du = -i(U, ζ* + f*) - d(p + k)
@Let dh = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))

#Checking the typings

#-----This part will be encapsulated into a function that automatically optimizes ----------
#Transforming the Forms expression into an Expression on arrays
exprs = [explicit(du); explicit(dh); explicit(pv)]
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
#println(rswstr)
#-------------------------------------------------------------------------------------------

#Testing the function

#Defining the Mesh
nx = 75
ny = 75
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

Lx, Ly = (1,1)
mesh = Arrays.Mesh(nx, ny, nh, msk, Lx, Ly)

#Initial Conditions
h0 = 0.15
H = 1
sigma = 0.05
gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))

h = zeros((nx,ny))

config = "dipole"

for i in nh+1:nx-nh, j in nh+1:ny-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]

	if config == "dipole"
		d=0.05

		h[i,j] = (H + h0 * (gaussian(x, y, 0.5+d/2, 0.5, sigma) - gaussian(x, y, 0.5-d/2, 0.5, sigma))) * mesh.A[5,5]
	elseif config == "vortex"
		h[i,j] = (H + h0 * gaussian(x, y, 0.5, 0.5, sigma)) * mesh.A[5,5]
	end
end

f =  15 .* ones((nx,ny)) .* mesh.A .* mesh.msk2d

u_x = zeros((nx, ny))
u_y = zeros((nx, ny))

zeta = zeros((nx,ny))
pv = zeros((nx,ny))

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
function rk3step!(dt, mesh, f, h, u_x, u_y, zeta, pv,
		dh1, dh2, dh3,
		du_x1, du_x2, du_x3,
		du_y1, du_y2, du_y3)
	
	rsw!(mesh ;f=f, h=h, u_x=u_x, u_y=u_y, zeta=zeta, pv=pv, dh=dh1, du_x=du_x1, du_y=du_y1)
	h .+= dt .* dh1
	u_x .+= dt .* du_x1
	u_y .+= dt .* du_y1

	rsw!(mesh ;f=f, h=h, u_x=u_x, u_y=u_y, zeta=zeta, pv=pv, dh=dh2, du_x=du_x2, du_y=du_y2)
	h .+= dt .* (-3/4 .* dh1 + 1/4 .* dh2)
	u_x .+= dt .* (-3/4 .* du_x1 + 1/4 .* du_x2)
	u_y .+= dt .* (-3/4 .* du_y1 + 1/4 .* du_y2)

	rsw!(mesh ;f=f, h=h, u_x=u_x, u_y=u_y, zeta=zeta, pv=pv, dh=dh3, du_x=du_x3, du_y=du_y3)
	h .+= dt .* (-1/12 .* dh1 - 1/12 .* dh2 + 2/3 .* dh3)
	u_x .+= dt .* (-1/12 .* du_x1 - 1/12 .* du_x2 + 2/3 .* du_x3)
	u_y .+= dt .* (-1/12 .* du_y1 - 1/12 .* du_y2 + 2/3 .* du_y3)
end

#Saving/Viz
save_every = 20

#Plotting
plot = true
var_to_plot = pv

if plot
	using Plots
	global hm = heatmap(var_to_plot[nh+2:nx-nh, nh+2:ny-nh])
	global anim = Animation()
	frame(anim, hm)
end

#NCDF
#TODO

#TimeLoop
tend = 1000
maxite = 10000
ite = 0

#todo "borrowed" from fluids2d, check further
cfl = 0.9
dtmax = 0.15
dt = dtmax
t = 0

tstart = time()
for ite in 1:maxite
	if t>=tend || dt<1e-5
		break
	end

	#Actual progress
	rk3step!(dt, mesh, f, h, u_x, u_y, zeta, pv,
		 dh1, dh2, dh3,
		 du_x1, du_x2, du_x3,
		 du_y1, du_y2, du_y3)

	global t+=dt
	maxU = maximum(abs.(u_x)/mesh.A[1,1])+maximum(abs.(u_y)/mesh.A[1,1])+1e-10
	global dt = min(cfl/maxU, dtmax)

	print("ite : $(ite)/$(maxite), dt: $(round(dt; digits = 2)), t : $(round(t; digits = 2))/$(tend)            \r")	
	
	if plot & ite%save_every==0
		global hm = heatmap(var_to_plot[nh+2:nx-nh, nh+2:ny-nh])
		frame(anim, hm)
	end
end

if plot
	mp4(anim, "out.mp4"; fps = 60)
end

println("")
println(time()-tstart)
