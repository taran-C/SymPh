include("../src/Poisson2D.jl")
using Plots

struct Mesh
	msk
	nx
	ny
	dx
	dy
end

nx = 100
ny = 100
nh = 3

msk = zeros(nx, ny)
msk[1+nh:nx-nh, 1+nh:ny-nh] .= 1


mesh = Mesh(msk, nx,ny,1/(nx-nh),1/(ny-nh))

poisson! = get_poisson_solver(mesh, "dirichlet")

x = zeros(nx, ny)
b = zeros(nx, ny)
b[nx÷2, ny÷3] = 1
b[nx÷2, 2*ny÷3] = -1

@time Main.poisson!(x,b)

hm = contour(x, fill=true)
display(hm)
