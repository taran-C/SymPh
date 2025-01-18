export Primal, Dual
export Form, FormVariable
export VectorVariable

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
	Form

D : degree (0/1/2... form)
P : primality
"""
abstract type Form{D,P<:Primality} end

struct FormVariable{D,P} <: Form{D,P}
	name::String
end

