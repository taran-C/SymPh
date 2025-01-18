# SymbolicPhysics
A module made to represent physical or mathematical operations on scalar/vector fields, and convert them into efficient kernels, to allow for lazy but powerful evaluation of operations.

## Structure
The module is composed of a few submodules, each corresponding to a different level of abstraction.
Here they are from most abstract to least abstract.

### Physics
Scalar, Vectors, associated with units to check consistency and classical physical operations (grad, div, curl, ...)
The `abstract` function allows us to convert into mathematical expressions.

### Maths
Vectors, Forms of different degrees and primalities, and their associated operations (exterior derivative, interior product, ...)
Using the `explicit` function returns us array expressions.

### Arrays


### Kernels
Functions that actually compute our operations on data.
