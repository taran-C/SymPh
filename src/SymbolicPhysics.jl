module SymbolicPhysics

import Base: +,*,^,-,/,<,>,!, string, getindex

#Arrays
module Arrays

include("Arrays/Variables.jl")
include("Arrays/Operators.jl")
include("Arrays/Conditionals.jl")

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
