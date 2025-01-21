# SymbolicPhysics
A module made to represent physical or mathematical operations on scalar/vector fields, and convert them into efficient kernels, to allow for lazy but powerful evaluation of operations.

## Structure
The module is composed of a few submodules, each corresponding to a different level of abstraction.
Here they are from most abstract to least abstract.

### Physics
Scalar, Vectors, associated with units to check consistency and classical physical operations (grad, div, curl, ...).

The spatial characteristics (metric, domain) should be in an ad-hoc object converted into maths.

The `abstract` function allows us to convert into mathematical expressions.

### Maths
Vectors, Forms of different degrees and primalities, and their associated operations (exterior derivative, interior product, ...).

Idem for space with a `Manifold` object.

Using the `explicit` function returns us array expressions.

### Arrays
Basic Algebraic expressions, which can be (thanks to `to_deptree` and `to_sequence`) encapsulated into `Block`s and `Sequence`s which respectively represent a number of equations to be computed simultaneously in our loop, or a sequence of blocks each to be computed in its own loop.

Here we directly use a `Mesh` object

Finally we use a `to_kernel` function to transform these expression on symbolic arrays into actual functions that use our data.

### Kernels
Functions that actually compute our operations on data.
