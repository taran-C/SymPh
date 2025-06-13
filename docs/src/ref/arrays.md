```@meta
CurrentModule = SymPh.Arrays
```

# Arrays
```@docs
Arrays
```

## Objects and Variables

### Abstract types
```@docs
Expression
Atom
Literal
Variable
```

### Concrete types
```@docs
RealValue
ScalarVariable
ArrayVariable
```

## Operations on arrays

```@docs
Operator
FuncCall
UnaryOperator
Negative
AbsoluteValue
BinaryOperator
Addition
Substraction
Multiplication
Division
Exponentiation
BooleanExpression
UnaryBooleanOperator
BinaryBooleanOperator
TernaryOperator
GreaterThan
```

## Meshes

```@docs
Mesh
CartesianMesh
PolarMesh
```

## Interpolations
```@docs
upwind
avg2pt
weno
```

## Finite volumes and finite differences
```@docs
fvtofd2
fvtofd4
fdtofv2
fdtofv4
```

## From arrays to kernels
```@docs
DepNode
Block
LoopBlock
CallBlock
Sequence
to_kernel
```
