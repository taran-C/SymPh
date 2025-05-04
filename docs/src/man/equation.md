# Defining our equation

## The maths
We will us as an example the stream function-vorticity formulation of the Euler equations with the vorticity $\omega$ as our prognostic variable. This gives us 
$ \partial_t \omega = \mathcal(L)_U \omega = \mathrm{d}(\iota_U \omega) $

with 

$ \begin{cases}
  \nabla^2 \omega^2 = \Psi^2 \\
  u^1 = \delta \Psi \\
  U = u^\sharp
  \end{cases}
 $

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

