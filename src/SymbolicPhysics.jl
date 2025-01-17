module SymbolicPhysics

import Base: +,*,^,-,/,<,>,!, string, getindex

#Arrays
include("Arrays/Variables.jl")
include("Arrays/Operators.jl")
include("Arrays/Conditionals.jl")

#Maths
include("Maths/Forms.jl")
include("Maths/DifferentialOperators.jl")

end
