# Defining our equation

## The maths
We will us as an example the stream function-vorticity formulation of the Euler equations with the vorticity $\omega$ as our prognostic variable. 

This gives us 

$$\partial_t \omega = \mathcal{L}_U \omega = \mathrm{d}(\iota_U \omega)$$

with 

$$\begin{cases}
    \nabla^2 \omega^2 = \Psi^2 \\
    u^1 = \delta \Psi \\
    U = u^\sharp
\end{cases}$$

## Representing it

We use the `@Let` macro to avoid having to specify the name of each object we define. It will use the variable name as the object name.

We start by defining our prognostic variable :

```julia
@Let omega = FormVariable{2, Dual}()
```

Then we define the rest of the variables in relation to already defined variables :

```julia
@Let psi = InverseLaplacian(omega)
@Let u = Codifferential(psi)
@Let U = Sharp(u)

#Time derivative
@Let dtomega = - ExteriorDerivative(InteriorProduct(U, omega)) #dtω = L(U,ω)
```

For now (if we want to use the time integrator from SymPh) the time derivative of prognostics variables must have the same name as the variable with `dt` at the beginning.

## Specifying the numerical methods

We then build an `ExplicitParam` object, specifying for example the interpolation method used for the interior product.

```julia
explparams = ExplicitParam(; interp = Arrays.weno)
```

## Building the kernel function

We can now use those objects to finally generate the function that will perform our computations.

```julia
euler_rhs! = to_kernel(dtomega; save = ["u_i", "u_j", "ι_U_omega_i", "ι_U_omega_j"], explparams = explparams)
```

We have to pass the objects that we want to be computed, for example in the case of a time integration, our time derivative, as well as a list of intermediary values that we want saved and our numerical methods.

## Building the mesh

We then define the parameters of a mesh, notably its size and a mask of its fluid cells

```julia
ni = 100
nj = 100
nh = 3

msk = zeros(ni, nj)
msk[nh+1:ni-nh, nh+1:nj-nh] .= 1

Lx, Ly = (1,1)
```

We also have to choose a loop manager to run the program

```julia
scalar = PlainCPU()
simd = VectorizedCPU(16)
threads = MultiThread(scalar)
```

And we can finally build our mesh object

```julia
mesh = Arrays.Mesh(ni, nj, nh, thsimd, msk, Lx, Ly)
```

## Defining the state and initial conditions

We initialize a `State` object, that will automatically allocate the arrays we ask it to give us, if they are not already initialized

```julia
state = State(mesh)
```

and we fill our $$\omega$$ with an initial condition, for example a vortex merging experiment :

```julia
gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
dipole(x,y,x0,y0,d,sigma) = gaussian(x, y, x0+d/2, y0, sigma) + gaussian(x, y, x0-d/2, y0, sigma)

omega = state.omega
for i in nh+1:nx-nh, j in nh+1:ny-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]
	omega[i,j] = tripole(x, y, 0.5,0.5,0.3,0.05) * mesh.msk2d[i,j]
end
```

We can now finally create the integration model that will hold all of this, by also specifying our prognostic variable and an integrator.

```julia
model = Model(euler_rhs!, mesh, state, ["omega"]; cfl = 100., dtmax = 5., integratorstep! = rk3step!)
```

## Running the simulation

We can now run our model for a certain time, and choose which variables should be written to disk.

```julia
run!(model; save_every = 5, profiling = false, tend = 10000, maxite = 100, writevars = (:u_i, :u_j, :omega, :psi))
```
