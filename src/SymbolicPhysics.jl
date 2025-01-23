module SymbolicPhysics

import Base: +,*,^,-,/,<,>,!, string, getindex

#Arrays
module Arrays
import Base: +,*,^,-,/,<,>,!, string, getindex

include("Arrays/Variables.jl")
include("Arrays/Operators.jl")
include("Arrays/DepTree.jl")
include("Arrays/ExecFlow.jl")
include("Arrays/ToKernels.jl")
end


#Maths
module Maths

import ..Arrays
import Base: +

include("Maths/Variables.jl")
include("Maths/Operators.jl")
include("Maths/Explicit.jl")
end

end
