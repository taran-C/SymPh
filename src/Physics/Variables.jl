abstract type Dim{T, L, M, I, Th, N, J} end

Time = Dim{1,0,0,0,0,0,0}
Mass = Dim{0,0,1,0,0,0,0}

Length = Dim{0,1,0,0,0,0,0}
Area = Dim{0,2,0,0,0,0,0}
Volume = Dim{0,3,0,0,0,0,0}

Speed = Dim{-1,1,0,0,0,0,0}
Acceleration = Dim{-2,1,0,0,0,0,0}
Force = Dim{-2,1,1,0,0,0,0}

mutable struct Variable{Dim}
	name::String
end

function unit(var::Variable{Dim{T, L, M, I, Th, N, J}}) where {T, L, M, I, Th, N, J}
	s = ""
	for (amt, name) in zip([T, L, M, I, Th, N, J],["T", "L", "M", "A", "ϴ", "N", "J"])
		if amt != 0
			s = s*"["*name*"]"*string(amt)
		end
	end
	println(s)
end
