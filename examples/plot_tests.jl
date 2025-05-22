using SymPh, SymPh.Maths
import SymPh.Arrays
using GLMakie

@Let u = FormVariable{1, Dual}()
@Let s = FormVariable{2, Primal}()

ni = 50
nj = 50
nh = 3

msk = zeros(ni, nj)
msk[nh+1:ni-nh, nh+1:nj-nh] .= 1

Lx = 1
Ly = 1

mesh = Arrays.PolarMesh(ni, nj, nh, Nothing, msk)
state = State(mesh)

xs = mesh.xc
ys = mesh.yc
state.u_i .= sin.(xs) .* cos.(ys) .* mesh.dx
state.u_j .= -cos.(xs) .* sin.(ys) .* mesh.dy
state.s .= sqrt.(state.u_j .^2 + state.u_i .^2)

inner = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)

f = Figure(size = (800, 800))
Axis(f[1, 1])#, backgroundcolor = "black")
plotform!(s, mesh, state)
plotform!(u, mesh, state)
display(f)
