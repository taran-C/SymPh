using SymbolicPhysics.Arrays

x = ArrayVariable("x")

y = x[0,1]+2*x[0,-1]
println(string(y))

f = to_kernels(y)

nx = 10
ny = 10
nh = 1

x = ones(nx,ny)
out = zeros(nx,ny)
display(x)
display(out)

f(;nx=nx, ny=ny, nh=nh, x=x, out=out)

display(x)
display(out)
