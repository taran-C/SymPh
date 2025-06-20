module SymPh

import Base: +,*,^,-,/,<,>,!, string, getindex

include("Misc.jl")
include("State.jl")
include("Integration.jl")
include("Model.jl")
include("plotform.jl")

"""
	Arrays

A module containing symbolic representations of arrays, and operations on them, and a way to convert them to computational kernels
"""
module Arrays
import Base: +,*,^,-,/,<,>,!, string, getindex, abs

include("Arrays/Variables.jl")
include("Arrays/Operators.jl")
include("Arrays/DepTree.jl")
include("Arrays/ExecFlow.jl")
include("Arrays/ToKernels.jl")
include("Arrays/Mesh.jl")
include("Arrays/Interpolations.jl")
include("Arrays/FVtoFD.jl")
end

"""
	Maths

A module containing symbolic representations of differential forms, and operations on them, and a way to convert them to array operations
"""
module Maths

import ..Arrays
import Base: +,-,*,/

include("Poisson2D.jl")
include("Maths/Variables.jl")
include("Maths/Operators.jl")
include("Maths/Explicit.jl")
include("Maths/ToKernel.jl")
end


end
