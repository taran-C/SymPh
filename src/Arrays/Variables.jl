export RealValue
export ScalarVariable, ArrayVariable

"""
	Expression

Generic Expression
"""
abstract type Expression end

getindex(A::Expression, depi, depj) = A
eval(expr::Expression) = eval(expr, Dict())

"""
	Atom <: Expression

An atom is a singular element holding a value (Variable or not)
"""
abstract type Atom <: Expression end
prec(expr::Atom) = 10

"""
	Literal <: Atom

Literal constant value
"""
abstract type Literal <: Atom end

"""
	RealValue(name::String, val::Real) <: Literal

A constant real value
"""
struct RealValue <: Literal
	name::String
	val::Real
end
string(expr::RealValue) = string(expr.val)
eval(expr::RealValue, vals::AbstractDict) = expr.val
RealValue(val::Real) = RealValue(string(val), val)

"""
	Variable <: Atom

Any variable value
"""
abstract type Variable <: Atom end

"""
	ScalarVariable(name::String) <: Variable

A variable holding a single scalar value
"""
struct ScalarVariable <: Variable
	name::String
end
string(expr::ScalarVariable) = expr.name
function eval(expr::ScalarVariable, vals::AbstractDict)
	if !haskey(vals, expr.name)
		return expr
	end
	return vals[expr.name]
end

"""
	ArrayVariable(name::String, depi::Integer, depj::Integer) <: Variable

Array object representing a variable name and a relative position
"""
struct ArrayVariable <: Variable
        name :: String
        depi :: Integer
        depj :: Integer
end
string(expr::ArrayVariable) = "$(expr.name)[$(expr.depi)+i,$(expr.depj)+j]"
ArrayVariable(name :: String) = ArrayVariable(name, 0, 0)
getindex(A::ArrayVariable, depi, depj) = ArrayVariable(A.name, A.depi+depi, A.depj+depj)
function eval(expr::ArrayVariable, vals::AbstractDict)
	if !haskey(vals, expr.name)
		return expr
	end
	if !haskey(vals, "i") | !haskey(vals, "j")
		throw(ErrorException("Can't get an array without coordinates"))
	end
	return vals[expr.name][vals["i"]+expr.depi, vals["j"]+expr.depj]
end



