export Primal, Dual
export Form, FormVariable

"""
	Primality

(Primal or Dual)
"""
abstract type Primality end
abstract type Primal <: Primality end
abstract type Dual <: Primality end

"""
	Form

D : degree (0/1/2... form)
P : primality
"""
abstract type Form{D,P} end

struct FormVariable{D,P} <: Form{D,P}
	name::String
end
