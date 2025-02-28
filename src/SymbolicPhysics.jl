module SymbolicPhysics

import Base: +,*,^,-,/,<,>,!, string, getindex

include("Misc.jl")

#Arrays
module Arrays
import Base: +,*,^,-,/,<,>,!, string, getindex

include("Arrays/Variables.jl")
include("Arrays/Operators.jl")
include("Arrays/DepTree.jl")
include("Arrays/ExecFlow.jl")
include("Arrays/ToKernels.jl")
include("Arrays/Mesh.jl")
include("Arrays/Interpolations.jl")
end


#Maths
module Maths

import ..Arrays
import Base: +,-,*,/

include("Maths/Variables.jl")
include("Maths/Operators.jl")
include("Maths/Explicit.jl")
end


end
