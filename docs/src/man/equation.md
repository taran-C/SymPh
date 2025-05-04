# Defining our equation

## The maths
We will us as an example the stream function-vorticity formulation of the Euler equations with the vorticity $\omega$ as our prognostic variable. This gives us $$\partial_t \omega = \mathcal(L)_U \omega = \mathrm{d}(\iota_U \omega)$$

with 

$$\begin{cases}
    \nabla^2 \omega^2 = \Psi^2 \\
    u^1 = \delta \Psi \\
    U = u^\sharp
\end{cases}$$

## Representing it

We use the `@Let` macro to avoid having to specify the name of each object we define. It will use the variable name as the object name.

We start by defining our prognostic variable :

```
@Let omega = FormVariable{2, Dual}()
```

Then we define the rest of the variables in relation to already defined variables :

```
@Let psi = InverseLaplacian(omega)
@Let u = Codifferential(psi)
@Let U = Sharp(u)

#Time derivative
@Let dtomega = - ExteriorDerivative(InteriorProduct(U, omega)) #dtω = L(U,ω)
```

For now (if we want to use the time integrator from SymPh) the time derivative of prognostics variables must have the same name as the variable with `dt` at the beginning.

## Specifying the numerical methods

We then build an `ExplicitParam` object, specifying for example the interpolation method used for the interior product.

```
explparams = ExplicitParam(; interp = Arrays.weno)
```

## Building the kernel function

We can now use those objects to finally generate the function that will perform our computations.

```
euler_rhs! = to_kernel(dtomega; save = ["u_x", "u_y", "ι_U_omega_x", "ι_U_omega_y"], explparams = explparams)
```

We have to pass the objects that we want to be computed, for example in the case of a time integration, our time derivative, as well as a list of intermediary values that we want saved and our numerical methods.
