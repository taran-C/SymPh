"""
Operators
"""

abstract type Operator <: Expression end

"""
Unary Operators
"""
abstract type UnaryOperator <: Operator end
getindex(expr::UnaryOperator, depx, depy) = op(expr)(expr.expr[depx,depy])
string(expr::UnaryOperator) = "$(symbol(expr))$(string(expr.expr))"
eval(expr::UnaryOperator, vals::AbstractDict) = op(expr)(eval(expr.expr, vals))

"""
Negative, represents the negation of an expression
"""
struct Negative{E<:Expression} <: UnaryOperator
	expr::E
end
symbol(expr::Negative) = "-"
op(expr::Negative) = -
prec(expr::Negative) = 2
-(expr::Expression) = Negative(expr)

"""
Binary Operators
"""
abstract type BinaryOperator <: Operator end
getindex(expr::BinaryOperator, depx, depy) = op(expr)(expr.left[depx,depy], expr.right[depx,depy])
function string(expr::BinaryOperator)
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
end
eval(expr::BinaryOperator, vals::AbstractDict) = op(expr)(eval(expr.left, vals), eval(expr.right, vals))

"""
Addition
"""
struct Addition{L<:Expression,R<:Expression} <: BinaryOperator
	left::L
	right::R
end
symbol(expr::Addition) = "+"
op(expr::Addition) = +
prec(expr::Addition) = 1
+(left::Expression, right::Expression) = Addition(left, right)
+(left::Expression, right::Real) = Addition(left, RealValue(right))
+(left::Real, right::Expression) = Addition(RealValue(left), right)

"""
Substraction
"""
struct Substraction{L<:Expression,R<:Expression} <: BinaryOperator
	left::L
	right::R
end
symbol(expr::Substraction) = "-"
op(expr::Substraction) = -
prec(expr::Substraction) = 2
-(left::Expression, right::Expression) = Substraction(left, right)
-(left::Expression, right::Real) = Substraction(left, RealValue(right))
-(left::Real, right::Expression) = Substraction(RealValue(left), right)

"""
Multiplication
"""
struct Multiplication{L<:Expression, R<:Expression} <: BinaryOperator
	left::L
	right::R
end
symbol(expr::Multiplication) = "*"
op(expr::Multiplication) = *
prec(expr::Multiplication) = 3
*(left::Expression, right::Expression) = Multiplication(left, right)
*(left::Expression, right::Real) = left * RealValue(right)
*(left::Real, right::Expression) = RealValue(left) * right

"""
Division

TODO error handling
"""
struct Division{L<:Expression, R<:Expression} <: BinaryOperator
	left::L
	right::R
end
symbol(expr::Division) = "/"
op(expr::Division) = /
prec(expr::Division) = 3
/(left::Expression, right::Expression) = Division(left, right)
/(left::Expression, right::Real) = Division(left, RealValue(right))
/(left::Real, right::Expression) = Division(RealValue(left), right)
