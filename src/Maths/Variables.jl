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
	VectorVariable
"""
abstract type Vect{P<:Primality} end

struct VectorVariable{P} <: Vect{P}
	name::String
end

"""
	Form{D,P<:Primality}

# Arguments
- `D::Integer` : degree (0/1/2... form)
- `P::Primality` : The primality of the form
"""
abstract type Form{D,P<:Primality} end

struct FormVariable{D,P} <: Form{D,P}
	name::String
end

degree(f::Form{D,P}) where {D,P} = D
primality(f::Form{D,P}) where {D,P} = P
primality(v::Vect{P}) where {P} = P
