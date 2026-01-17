export Primal, Dual
export Form, FormVariable
export VectorVariable
export Vect
export degree, primality

"""
	Primality

(Primal or Dual)
"""
abstract type Primality end
abstract type Primal <: Primality end
abstract type Dual <: Primality end

"""
	Vect{P<:Primality}

A vector field
"""
abstract type Vect{P<:Primality} end

"""
	VectorVariable{P}(name::String, save::Bool) <: Vect{P}

A variable representing a named vector field
"""
struct VectorVariable{P} <: Vect{P}
	name::String
	save::Bool
end
VectorVariable{P}(; name = "", save=false) where {P} = VectorVariable{P}(name, save)

"""
	Form{D,P<:Primality}

# Arguments
- `D::Integer` : degree (0/1/2... form)
- `P::Primality` : The primality of the form
"""
abstract type Form{D,P<:Primality} end

"""
	FormVariable{D,P}(name::String, save::Bool) <: Form{D,P}

A variable representing a named differential form
"""
struct FormVariable{D,P} <: Form{D,P}
	name::String
	save::Bool
end
FormVariable{D,P}(; name = "", save=false) where {D,P} = FormVariable{D,P}(name, save)

degree(f::Form{D,P}) where {D,P} = D
primality(f::Form{D,P}) where {D,P} = P
primality(v::Vect{P}) where {P} = P
