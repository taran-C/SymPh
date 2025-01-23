export RealValue
export ScalarVariable, ArrayVariable

"""
Generic Expression
"""
abstract type Expression end

getindex(A::Expression, depx, depy) = A
eval(expr::Expression) = eval(expr, Dict())

"""
An atom is a singular element holding a value (Variable or not)
"""
abstract type Atom <: Expression end
prec(expr::Atom) = 10

"""
Literals (literal constant value)
"""
abstract type Literal <: Atom end

struct RealValue <: Literal
	name::String
	val::Real
end
string(expr::RealValue) = string(expr.val)
eval(expr::RealValue, vals::AbstractDict) = expr.val
RealValue(val::Real) = RealValue(string(val), val)

"""
Variables
"""
abstract type Variable <: Atom end

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

#Array object representing a variable name and a relative position
struct ArrayVariable <: Variable
        name :: String
        depx :: Integer
        depy :: Integer
end
string(expr::ArrayVariable) = "$(expr.name)[$(expr.depx)+i,$(expr.depy)+j]"
ArrayVariable(name :: String) = ArrayVariable(name, 0, 0)
getindex(A::ArrayVariable, depx, depy) = ArrayVariable(A.name, A.depx+depx, A.depy+depy)
function eval(expr::ArrayVariable, vals::AbstractDict)
	if !haskey(vals, expr.name)
		return expr
	end
	if !haskey(vals, "i") | !haskey(vals, "j")
		throw(ErrorException("Can't get an array without coordinates"))
	end
	return vals[expr.name][vals["i"]+expr.depx, vals["j"]+expr.depy]
end



