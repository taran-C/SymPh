```@meta
CurrentModule = SymPh.Maths
```
# Differential Geometry

## Abstract Objects

These are abstract object representing differential forms and vectors

```@docs
Form
Vect
```

## Variables

The variables upon which we act with our operators. Typically our prognostic variables will be defined like this since they are not defined as an operation on something.

```@docs
FormVariable
VectorVariable
```

## Operators

These are the operators acting on our variables. Their type is the type of the result of the operation.

### The Wedge product $$\wedge$$
```@docs
Wedge
```

### The Exterior Derivative $$\mathrm{d}$$
```@docs
ExteriorDerivative
```

The exterior derivative $$\mathrm{d}\omega$$ of a $$k$$-form $$\omega$$ is a $$k+1$$-form.

In 2D, 

$$\mathrm{d}F^0 = \partial_x f \mathrm{d} x+\partial_y f \mathrm{d} y$$

Corresponds to the gradient operation on a scalar field.

$$\mathrm{d}U^1 = (\partial_x u_2-\partial_y u_1) \mathrm{d} x \wedge \mathrm{d} y$$

Corresponds to the curl of a vector field.

$$\mathrm{d} W^2 = 0$$

### The Codifferential $$\delta$$
```@docs
Codifferential
```

### The Interior Product with a vector field $$\iota_\mathbf{X}$$
```@docs
InteriorProduct
```

Let $$\mathbf{X}$$ be a vector field. Then :

$\iota_\mathbf{X} : \Omega^k (M) \mapsto \Omega^{k-1} (M)$

We have $$\iota_\mathbf{X} (\alpha)= \alpha (\mathbf{X})$$, and $$\iota_\mathbf{X} ( \alpha \wedge \beta) = ( \iota_\mathbf{X}  \alpha) \wedge  \beta + (-1)^{|\alpha|}  \alpha \wedge ( \iota_\mathbf{X} \beta )$$ with $$|\alpha|=k$$ the degree of the $$k$$-form $$\alpha$$.


Let $\mathbf{X}=x_1 \mathbf{i} + x_2 \mathbf{j}$, then (in 2D),

$$\iota_\mathbf{X} F^0 = 0$$

$$\iota_\mathbf{X} U^1 = \iota_\mathbf{X}(u_1 \mathrm{d} x+u_2 \mathrm{d} y)=u_1 \iota_\mathbf{X} \mathrm{d} x+u_2 \iota_\mathbf{X} \mathrm{d} y = u_1x_1+u_2 x_2$$

Corresponds to scalar product of two vector fields.

$$\iota_\mathbf{X} W^2 = \iota_\mathbf{X}(w \mathrm{d} x \wedge \mathrm{d} y) = w( \iota_\mathbf{X}( \mathrm{d} x) \wedge \mathrm{d} y - \mathrm{d} x \wedge \iota_\mathbf{X}( \mathrm{d}  y)) = w(x_1 \mathrm{d}  y-x_2 \mathrm{d}  x)$$

Corresponds to a kind of "perpendicular product" of a scalar field with a vector field.

### The Sharp $$\sharp$$
```@docs
Sharp
```

### The Hodge star $$\star$$
```@docs
Hodge
```

### The inverse of a Laplacian
```@docs
InverseLaplacian
```

### Others

```@docs
Addition
Substraction
Negative
FuncCall
Division
RealProdForm
```
