export Addition
export Substraction
export ExteriorDerivative
export InteriorProduct
export Sharp
export Hodge
export InnerProduct
export FuncCall

"""
	FuncCall
"""
mutable struct FuncCall{D, P} <: Form{D,P}
	name::String
	func
	args::Vector{Form}
end

"""
	Addition
"""
mutable struct Addition{D,P} <: Form{D,P}
	name::String
	left::Form{D,P}
	right::Form{D,P}
end
+(name::String, left::Form, right::Form) = Addition(name, left, right)
+(left::Form, right::Form) = Addition("P_"*left.name*"_"*right.name, left, right)

"""
	Substraction
"""
mutable struct Substraction{D,P} <: Form{D,P}
	name::String
	left::Form{D,P}
	right::Form{D,P}
end
-(name::String, left::Form, right::Form) = Substraction(name, left, right)
-(left::Form{D,P}, right::Form{D,P}) where {D,P} = Substraction{D,P}("P_"*left.name*"_"*right.name, left, right)

"""
	Negative
"""
mutable struct Negative{D,P} <: Form{D,P}
	name::String
	form::Form{D,P}
end
-(name::String, form::Form) = Negative(name, form)
-(form::Form) = Negative("N_"*form.name, form)

"""
	Division
	TODO What is that actually in terms of FORMS ?
	Simple division by a scalar field proxied by a 0-form ?
"""
mutable struct Division{D,P} <: Form{D,P}
	name::String
	left::Form{D,P}
	right::Form
end
/(name::String, left::Form, right::Form) = Division(name, left, right)


"""
	ExteriorDerivative
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

"""
	Sharp
"""
mutable struct Sharp{P} <: Vect{P}
	name::String
	form::Form{1, P}
end
Sharp(form::Form) = Sharp("#_"*form.name, form)

"""
	Hodge
"""
mutable struct Hodge{D, P} <: Form{D, P}
	name::String
	form::Form

	function Hodge(name::String, form::Form{D, P}) where {D, P}
		if P<:Primal
			return new{2-D, Dual}(name, form)
		else
			return new{2-D, Primal}(name, form)
		end
	end
	Hodge(form) = Hodge("*_"*form.name, form)
end

"""
	InnerProduct
"""
mutable struct InnerProduct{D, P} <: Form{D, P}
	name::String
	left::Form
	right::Form

	function InnerProduct(name::String, left::Form{D1, P}, right::Form{D2, P}) where {D1, D2, P}
		if P<:Primal
			return new{D1+D2, Dual}(name, left, right)
		else	
			return new{D1+D2, Primal}(name, left, right)
		end
	end
	InnerProduct(left, right) = InnerProduct("IP_"*left.name*"_"*right.name, left, right)
end

"""
	RealProducts
"""
mutable struct RealProdForm{D, P} <: Form{D, P}
	name::String
	real::Real
	form::Form{D, P}
end
*(name::String, real::Real, form::Form) = RealProdForm(name, real, form)
*(name::String, form::Form, real::Real) = RealProdForm(name, real, form)
*(real::Real,  form::Form) = RealProdForm("T_"*string(real)*"_"*form.name, real, form)
*(form::Form, real::Real) = RealProdForm("T_"*string(real)*"_"*form.name, real, form)
