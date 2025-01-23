"""
Operators
"""
abstract type Operator <: Expression end

"""
Unary Operators
"""
abstract type UnaryOperator <: Operator end
getindex(expr::UnaryOperator, depx, depy) = op(expr)(expr.expr[depx,depy])
string(expr::UnaryOperator; fpar=false) = fpar ? "($(symbol(expr))$(string(expr.expr)))" : "$(symbol(expr))$(string(expr.expr))"
eval(expr::UnaryOperator, vals::AbstractDict) = op(expr)(eval(expr.expr, vals))

"""
Negative, represents the negation of an expression
"""
struct Negative <: UnaryOperator
	name::String
	expr::Expression
end
symbol(expr::Negative) = "-"
op(expr::Negative) = -
prec(expr::Negative) = 2
-(expr::Expression) = Negative(expr.name*"_neg", expr)

"""
Binary Operators
"""
abstract type BinaryOperator <: Operator end
getindex(expr::BinaryOperator, depx, depy) = op(expr)(expr.left[depx,depy], expr.right[depx,depy])
function string(expr::BinaryOperator; fpar=false)
	if fpar==false
		ret = ""

		if prec(expr.left)<prec(expr)
			ret = string(ret, "($(string(expr.left)))")
		else
			ret = string(ret, string(expr.left))
		end
		
		ret = string(ret, " $(symbol(expr)) ")
		
		if prec(expr.right)<prec(expr)
			ret = string(ret, "($(string(expr.right)))")
		else
			ret = string(ret, string(expr.right))
		end
		
		return ret
	else
		return "($(string(expr.left)))$(symbol(expr))($(string(expr.right)))"
	end	
end
eval(expr::BinaryOperator, vals::AbstractDict) = op(expr)(eval(expr.left, vals), eval(expr.right, vals))

"""
Addition
"""
struct Addition <: BinaryOperator
	name::String
	left::Expression
	right::Expression
end
symbol(expr::Addition) = "+"
op(expr::Addition) = +
prec(expr::Addition) = 1
+(left::Expression, right::Expression) = Addition("p_"*left.name*"_"*right.name, left, right)
+(left::Expression, right::Real) = Addition(left, RealValue(right))
+(left::Real, right::Expression) = Addition(RealValue(left), right)

"""
Substraction
"""
struct Substraction <: BinaryOperator
	name::String
	left::Expression
	right::Expression
end
symbol(expr::Substraction) = "-"
op(expr::Substraction) = -
prec(expr::Substraction) = 2
-(left::Expression, right::Expression) = Substraction("m_"*left.name*"_"*right.name, left, right)
-(left::Expression, right::Real) = Substraction(left, RealValue(right))
-(left::Real, right::Expression) = Substraction(RealValue(left), right)

"""
Multiplication
"""
struct Multiplication <: BinaryOperator
	name::String
	left::Expression
	right::Expression
end
symbol(expr::Multiplication) = "*"
op(expr::Multiplication) = *
prec(expr::Multiplication) = 3
*(left::Expression, right::Expression) = Multiplication("t_"*left.name*"_"*right.name, left, right)
*(left::Expression, right::Real) = left * RealValue(right)
*(left::Real, right::Expression) = RealValue(left) * right

"""
Division

TODO error handling
"""
struct Division <: BinaryOperator
	name::String
	left::Expression
	right::Expression
end
symbol(expr::Division) = "/"
op(expr::Division) = /
prec(expr::Division) = 3
/(left::Expression, right::Expression) = Division("d_"*left.name*"_"*right.name, left, right)
/(left::Expression, right::Real) = Division(left, RealValue(right))
/(left::Real, right::Expression) = Division(RealValue(left), right)


abstract type BinaryBooleanOperator <: BinaryOperator end
abstract type UnaryBooleanOperator <: UnaryOperator end

BooleanExpression = Union{UnaryBooleanOperator, BinaryBooleanOperator}

export TernaryOperator
"""
TernaryOperator

	symbolic representation of a ? b : c
"""
struct TernaryOperator <: Operator
	name::String
	a::BooleanExpression
	b::Expression
	c::Expression
end
#TODO check with true if else cause i can't seem to find a way to override the ternary operator ?
eval(expr::TernaryOperator, vals::AbstractDict) = eval(expr.a, vals) ? eval(expr.b, vals) : eval(expr.c, vals)
string(expr::TernaryOperator) = "($(string(expr.a))) ? ($(string(expr.b))) : ($(string(expr.c)))"

#Conditionals
"""
GreaterThan

	tests if left > right
"""
struct GreaterThan <: BinaryBooleanOperator
	name::String
	left::Expression
	right::Expression
end
symbol(expr::GreaterThan) = ">"
op(expr::GreaterThan) = >
prec(expr::GreaterThan) = 10
>(left::Expression, right::Expression) = GreaterThan(left.name*"_"*right.name*"_gt", left, right)
>(left::Expression, right::Real) = GreaterThan(left, RealValue(right))
>(left::Real, right::Expression) = GreaterThan(RealValue(left), right)

