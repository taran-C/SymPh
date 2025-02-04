export Addition
export ExteriorDerivative
export InteriorProduct

"""
	Addition
"""
mutable struct Addition{D,P} <: Form{D,P}
	name::String
	left::Form{D,P}
	right::Form{D,P}
end
+(left::Form{D,P}, right::Form{D,P}) where {D,P} = Addition{D,P}("P_"*left.name*"_"*right.name, left, right)

"""
	InteriorProduct
"""
mutable struct ExteriorDerivative{D,P} <: Form{D,P}
	name::String
	form::Form
	
	function ExteriorDerivative(name::String, expr::Form{D,P}) where {D,P}
		return new{D+1, P}(name, expr)
	end
	ExteriorDerivative(expr::Form) = ExteriorDerivative("d"*expr.name, expr)
end

"""
	InteriorProduct

TODO check primality (should it actually be implemented here or into Arrays ?
"""
mutable struct InteriorProduct{D, Pv, Pf} <: Form{D,Pf}
	name::String
	vect::Vect
	form::Form
	
	function InteriorProduct(name::String, vect::Vect{Pv}, form::Form{D,Pf}) where {Pv,D,Pf}
		return new{D-1, Pv, Pf}(name, vect, form)
	end
	InteriorProduct(vect::Vect, form::Form) = InteriorProduct("ι_"*vect.name*"_"*form.name, vect, form)
end
