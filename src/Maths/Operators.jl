export Addition
export Substraction
export ExteriorDerivative
export InteriorProduct
export Sharp
export Hodge
export InnerProduct
export FuncCall
export Codifferential
export InverseLaplacian
#TODO
#export InverseHodge
#export Flat

"""
	FuncCall{D, P}(name::String, func, args::Vector{Form}) <: Form{D,P}

Represents a call to `func` applied to the objects represented by `args`. Forces the computation of its arguments.
"""
mutable struct FuncCall{D, P} <: Form{D,P}
	name::String
	func
	args::Vector{Form}
end

"""
	Addition{D,P}(name::String, left::Form{D,P}, right::Form{D,P}) <: Form{D,P}

The element-wise sum ``left + right``.
"""
mutable struct Addition{D,P} <: Form{D,P}
	name::String
	left::Form{D,P}
	right::Form{D,P}
end
+(name::String, left::Form, right::Form) = Addition(name, left, right)
+(left::Form, right::Form) = Addition("P_"*left.name*"_"*right.name, left, right)

"""
	Substraction{D,P}(name::String, left::Form{D,P}, right::Form{D,P}) <: Form{D,P}

The substraction ``left - right``.
"""
mutable struct Substraction{D,P} <: Form{D,P}
	name::String
	left::Form{D,P}
	right::Form{D,P}
end
-(name::String, left::Form, right::Form) = Substraction(name, left, right)
-(left::Form{D,P}, right::Form{D,P}) where {D,P} = Substraction{D,P}("P_"*left.name*"_"*right.name, left, right)

"""
	Negative{D,P}(name::String, form::Form{D,P}) <: Form{D,P}

The inverse of a form.
"""
mutable struct Negative{D,P} <: Form{D,P}
	name::String
	form::Form{D,P}
end
-(name::String, form::Form) = Negative(name, form)
-(form::Form) = Negative("N_"*form.name, form)

"""
	Division{D,P}(name::String, left::Form{D,P}, right::Form) <: Form{D,P}
	
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
	ExteriorDerivative{D,P}(name::String, omega::Form{D-1,P}) <: Form{D,P}

The exterior derivative ``\\mathrm{d}\\omega``

Transforms a ``k``-form into a ``k+1``-form
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
	Codifferential{D,P}(name::String, omega::Form{D+1,P}) <: Form{D,P}

The codifferential ``\\delta \\omega``

Transforms a ``k``-form into a ``k-1``-form
"""
mutable struct Codifferential{D,P} <: Form{D,P}
	name::String
	form::Form

	function Codifferential(name::String, expr::Form{D,P}) where {D,P}
		return new{D-1, P}(name, expr)
	end
	Codifferential(expr::Form) = Codifferential("CODIF_"*expr.name, expr)
end


"""
	InteriorProduct{D, Pv, Pf}(name::String, X::Vect, omega::Form, interp = Nothing) <: Form{D, Pf}

The contraction of a ``k``-form ``\\omega`` with a vector field ``\\mathbf{X}`` which gives us a ``k-1``-form ``\\iota_\\mathbf{X}\\omega``

Possibility to specify interpolation function
"""
mutable struct InteriorProduct{D, Pv, Pf} <: Form{D,Pf}
	name::String
	vect::Vect
	form::Form
	interp
	
	function InteriorProduct(name::String, vect::Vect{Pv}, form::Form{D,Pf}; interp = Nothing) where {Pv,D,Pf}
		return new{D-1, Pv, Pf}(name, vect, form, interp)
	end
	InteriorProduct(vect::Vect, form::Form; interp = Nothing) = InteriorProduct("ι_"*vect.name*"_"*form.name, vect, form; interp=interp)
end

"""
	Sharp{P}(name::String, form::Form{1, P}) <: Vect{P}

Corresponds to an application of the metric
"""
mutable struct Sharp{P} <: Vect{P}
	name::String
	form::Form{1, P}
end
Sharp(form::Form) = Sharp("#_"*form.name, form)

"""
	Hodge{D, P}(name::String, form::Form) <: Form{D, P}

Brings a  ``k``-form to a ``n-k`` form and goes from dual to primal and inversely
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

#= TODO DO WE NEED IT ? DON'T THINK SO
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
=#

"""
	RealProdForm{D, P}(name::String, real::Real, form::Form) <: Form{D, P}

`real` time `form`
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

"""
	InverseLaplacian

Represents the solution to a Poisson problem
"""
mutable struct InverseLaplacian{D,P} <: Form{D,P}
	name::String
	form::Form{D,P}
end
InverseLaplacian(form::Form) = InverseLaplacian("INVLAP_"*form.name, form)
