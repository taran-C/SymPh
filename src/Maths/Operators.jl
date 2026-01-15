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
export Wedge
#TODO
#export InverseHodge
#export Flat

"""
	FuncCall{D, P}(name::String, func, args::Vector{Form}) <: Form{D,P}

Represents a call to `func` applied to the objects represented by `args`. Forces the computation of its arguments.
"""
mutable struct FuncCall{D, P} <: Form{D,P}
	name::String
	save::Bool
	func
	args::Vector{Form}
end
FuncCall{D,P}(func, args::Vector{Form}; name::String, save::Bool) where {D,P} = FuncCall{D,P}(name, save, func, args)

"""
	Addition{D,P}(name::String, left::Form{D,P}, right::Form{D,P}) <: Form{D,P}

The element-wise sum ``left + right``.
"""
mutable struct Addition{D,P} <: Form{D,P}
	name::String
	save::Bool
	left::Form{D,P}
	right::Form{D,P}
end
+(left::Form{D,P}, right::Form{D,P}; name="P_"*left.name*"_"*right.name, save=false) where {D,P} = Addition{D,P}(name, save, left, right)

"""
	Substraction{D,P}(name::String, left::Form{D,P}, right::Form{D,P}) <: Form{D,P}

The substraction ``left - right``.
"""
mutable struct Substraction{D,P} <: Form{D,P}
	name::String
	save::Bool
	left::Form{D,P}
	right::Form{D,P}
end
-(left::Form{D,P}, right::Form{D,P}; name="M_"*left.name*"_"*right.name, save=false) where {D,P} = Substraction{D,P}(name, save, left, right)

"""
	Negative{D,P}(name::String, form::Form{D,P}) <: Form{D,P}

The inverse of a form.
"""
mutable struct Negative{D,P} <: Form{D,P}
	name::String
	save::Bool
	form::Form{D,P}
end
-(form::Form; name="N_"*form.name, save=false) = Negative(name, save, form)

"""
	Division{D,P}(name::String, left::Form{D,P}, right::Form) <: Form{D,P}
	
TODO What is that actually in terms of FORMS ?
Simple division by a scalar field proxied by a 0-form ?
"""
mutable struct Division{D,P} <: Form{D,P}
	name::String
	save::Bool
	left::Form{D,P}
	right::Form
end
/(left::Form, right::Form; name="DIV_"*left.name*"_"*right.name, save=false) = Division(name, save, left, right)

"""
	Wedge{Dl + Dr,P}(name::String, left::Form{Dl,P}, right::Form{Dr,P}) <: Form{Dl + Dr, P}

Wedge product of two forms TODO Wedge between different primalities ?
"""
mutable struct Wedge{D, Dl, Dr, P} <: Form{D, P}
	name::String
	save::Bool
	left::Form
	right::Form

	function Wedge(name::String, save::Bool, left::Form{Dl, P}, right::Form{Dr, P}) where {Dl, Dr, P}
		return new{Dl+Dr, Dl, Dr, P}(name, save, left, right)
	end
end
Wedge(left::Form, right::Form; name="WEDGE_"*left.name*"_"*right.name, save=false) = Wedge(name, save, left, right)

"""
	ExteriorDerivative{D,P}(name::String, omega::Form{D-1,P}) <: Form{D,P}

The exterior derivative ``\\mathrm{d}\\omega``

Transforms a ``k``-form into a ``k+1``-form
"""
mutable struct ExteriorDerivative{D,P} <: Form{D,P}
	name::String
	save::Bool
	form::Form
	
	function ExteriorDerivative(name::String, save::Bool, expr::Form{D,P}) where {D,P}
		return new{D+1, P}(name, save, expr)
	end
end
ExteriorDerivative(expr::Form; name="d"*expr.name, save=false) = ExteriorDerivative(name, save, expr)

"""
	Codifferential{D,P}(name::String, omega::Form{D+1,P}) <: Form{D,P}

The codifferential ``\\delta \\omega``

Transforms a ``k``-form into a ``k-1``-form
"""
mutable struct Codifferential{D,P} <: Form{D,P}
	name::String
	save::Bool
	form::Form

	function Codifferential(name::String, save::Bool, expr::Form{D,P}) where {D,P}
		return new{D-1, P}(name, save, expr)
	end
end
Codifferential(expr::Form; name="CODIF_"*expr.name, save=false) = Codifferential(name, save, expr)

"""
	InteriorProduct{D, Pv, Pf}(name::String, X::Vect, omega::Form, interp = Nothing) <: Form{D, Pf}

The contraction of a ``k``-form ``\\omega`` with a vector field ``\\mathbf{X}`` which gives us a ``k-1``-form ``\\iota_\\mathbf{X}\\omega``

Possibility to specify interpolation function
"""
mutable struct InteriorProduct{D, Pv, Pf} <: Form{D,Pf}
	name::String
	save::Bool
	vect::Vect
	form::Form
	interp
	
	function InteriorProduct(name::String, save::Bool, vect::Vect{Pv}, form::Form{D,Pf}, interp) where {Pv,D,Pf}
		return new{D-1, Pv, Pf}(name, save, vect, form, interp)
	end
end
InteriorProduct(vect::Vect, form::Form; interp = nothing, name="ι_"*vect.name*"_"*form.name, save=false) = InteriorProduct(name, save, vect, form, interp)

"""
	Sharp{P}(name::String, form::Form{1, P}) <: Vect{P}

Corresponds to an application of the metric
"""
mutable struct Sharp{P} <: Vect{P}
	name::String
	save::Bool
	form::Form{1, P}
	fvtofd
	fdtofv
end
Sharp(form::Form; name="#_"*form.name, save=false, fvtofd=nothing, fdtofv=nothing) = Sharp(name, save, form, fvtofd, fdtofv)

"""
	Hodge{D, P}(name::String, form::Form) <: Form{D, P}

Brings a  ``k``-form to a ``n-k`` form and goes from dual to primal and inversely
"""
mutable struct Hodge{D, P} <: Form{D, P}
	name::String
	save::Bool
	form::Form

	function Hodge(name::String, save::Bool, form::Form{D, P}) where {D, P}
		if P<:Primal
			return new{2-D, Dual}(name, save, form)
		else
			return new{2-D, Primal}(name, save, form)
		end
	end
end
Hodge(form::Form; name="*_"*form.name, save=false) = Hodge(name, save, form)

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
	save::Bool
	real::Real
	form::Form{D, P}
end
*(real::Real, form::Form; name="T_"*string(real)*"_"*form.name, save=false) = RealProdForm(name, save, real, form)
*(form::Form, real::Real; name="T_"*string(real)*"_"*form.name, save=false) = RealProdForm(name, save, real, form)

"""
	InverseLaplacian

Represents the solution to a Poisson problem
"""
mutable struct InverseLaplacian{D,P} <: Form{D,P}
	name::String
	save::Bool
	form::Form{D,P}
end
InverseLaplacian(form::Form; name="INVLAP_"*form.name, save=true) = InverseLaplacian(name, save, form)
