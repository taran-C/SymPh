"""
	Operator <: Expression
	
Generic operation on an expression
"""
abstract type Operator <: Expression end

export FuncCall
"""
	FuncCall(name::String, func, args::Vector{Expression}, depi::Integer, depj::Integer) <: Operator

Forces computation of its arguments, then calls a function on them
For now the function must have signature `f(out, args...)` where out will be replaced by the name assigned to the operator
"""
mutable struct FuncCall <: Operator
	name::String
	save::Bool
	func#Figure out typing
	args::Vector{Expression}
	depi::Integer
	depj::Integer
end
FuncCall(name, func, args) = FuncCall(name, true, func, args, 0, 0)
getindex(expr::FuncCall, depi, depj) = FuncCall(expr.name, expr.save, expr.func, expr.args, expr.depi + depi, expr.depj + depj)
function string(expr::FuncCall)
	str = expr.name * "("
	for arg in expr.args
		str = str * string(arg) * ", "
	end
	str = chop(str, head = 0, tail = 2)

	return str
end

"""
	UnaryOperator
	
Any operator of a single argument
"""
abstract type UnaryOperator <: Operator end
getindex(expr::UnaryOperator, depi, depj) = typeof(expr)(expr.name, expr.save, expr.expr, expr.depi+depi, expr.depj+depj)
string(expr::UnaryOperator) = "($(symbol(expr))($(string(expr.expr[expr.depi, expr.depj]))))"
eval_expr(expr::UnaryOperator, vals::AbstractDict) = op(expr)(eval_expr(expr.expr, vals))

"""
	Negative
	
Represents the negation of an expression
"""
mutable struct Negative <: UnaryOperator
	name::String
	save::Bool
	expr::Expression
        depi::Integer
        depj::Integer
end
symbol(expr::Negative) = "-"
op(expr::Negative) = -
prec(expr::Negative) = 2
Negative(expr::Expression; name="neg_"*expr.name, save=false) = Negative(name, save, expr, 0, 0)
-(expr::Expression) = Negative(expr)

"""
	AbsoluteValue
"""
mutable struct AbsoluteValue <: UnaryOperator
	name::String
	save::Bool
	expr::Expression
	depi::Integer
	depj::Integer
end
symbol(expr::AbsoluteValue) = "abs"
op(expr::AbsoluteValue) = abs
prec(expr::AbsoluteValue) = 10
AbsoluteValue(expr::Expression; name="abs_"*expr.name, save=false) = AbsoluteValue(name, save, expr, 0, 0)
abs(expr::Expression) = AbsoluteValue(expr)

"""
	BinaryOperator

Any operator on two expressions
"""
abstract type BinaryOperator <: Operator end
getindex(expr::BinaryOperator, depi, depj) = typeof(expr)(expr.name, expr.save, expr.left, expr.right, expr.depi+depi, expr.depj+depj)
string(expr::BinaryOperator) = "($(string(expr.left[expr.depi, expr.depj])))$(symbol(expr))($(string(expr.right[expr.depi, expr.depj])))"
eval_expr(expr::BinaryOperator, vals::AbstractDict) = op(expr)(eval_expr(expr.left, vals), eval_expr(expr.right, vals))

"""
	Addition
"""
mutable struct Addition <: BinaryOperator
	name::String
	save::Bool
	left::Expression
	right::Expression
        depi::Integer
        depj::Integer
end
symbol(expr::Addition) = "+"
op(expr::Addition) = +
prec(expr::Addition) = 1
Addition(left::Expression, right::Expression; name="p_"*string(left)*"_"*right.name, save=false) = Addition(name, save, left, right, 0,0)
+(left::Expression, right::Expression) = Addition(left, right)
+(left::Real, right::Expression) = Addition(left, right)
+(left::Expression, right::Real) = Addition(left, right)

"""
	Substraction
"""
mutable struct Substraction <: BinaryOperator
	name::String
	save::Bool
	left::Expression
	right::Expression
        depi :: Integer
        depj :: Integer
end
symbol(expr::Substraction) = "-"
op(expr::Substraction) = -
prec(expr::Substraction) = 2
Substraction(left::Expression, right::Expression; name="m_"*left.name*"_"*right.name, save=false) = Substraction(name, save, left, right, 0,0)
-(left::Expression, right::Expression) = Substraction(left, right)
-(left::Expression, right::Real) = Substraction(left, RealValue(right))
-(left::Real, right::Expression) = Substraction(RealValue(left), right)

"""
	Multiplication
"""
mutable struct Multiplication <: BinaryOperator
	name::String
	save::Bool
	left::Expression
	right::Expression
        depi :: Integer
        depj :: Integer
end
symbol(expr::Multiplication) = "*"
op(expr::Multiplication) = *
prec(expr::Multiplication) = 3
Multiplication(left::Expression, right::Expression; name="t_"*left.name*"_"*right.name, save=false) = Multiplication(name, save, left, right, 0,0)
*(left::Expression, right::Expression) = Multiplication(left, right)
*(left::Real, right::Expression) = Multiplication(RealValue(left), right)
*(left::Expression, right::Real) = Multiplication(left, RealValue(right))

"""
	Division

TODO error handling
"""
mutable struct Division <: BinaryOperator
	name::String
	save::Bool
	left::Expression
	right::Expression
        depi :: Integer
        depj :: Integer
end
symbol(expr::Division) = "/"
op(expr::Division) = /
prec(expr::Division) = 3
Division(left::Expression, right::Expression; name = "d_"*left.name*"_"*right.name, save = false) = Division(name, save, left, right, 0,0)
/(left::Expression, right::Expression) = Division(left, right)
/(left::Expression, right::Real) = left / RealValue(right)
/(left::Real, right::Expression) = RealValue(left) / right

"""
	Exponentiation
"""
mutable struct Exponentiation <: BinaryOperator
	name::String
	save::Bool
	left::Expression
	right::Expression
	depi::Integer
	depj::Integer
end
symbol(expr::Exponentiation) = "^"
op(expr::Exponentiation) = ^
prec(expr::Exponentiation) = 10
Exponentiation(left, right; name="pow_"*left.name*"_"*right.name, save=false) = Exponentiation(name, save, left, right, 0,0)
^(left::Expression, right::Expression) = Exponentiation(left, right)
^(left::Expression, right::Real) = left ^ RealValue(right)
^(left::Real, right::Expression) = RealValue(left) ^ right

"""
	BinaryBooleanOperator
"""
abstract type BinaryBooleanOperator <: BinaryOperator end
"""
	UnaryBooleanOperator
"""
abstract type UnaryBooleanOperator <: UnaryOperator end

"""
	BooleanExpression
"""
BooleanExpression = Union{UnaryBooleanOperator, BinaryBooleanOperator}

export TernaryOperator
"""
	TernaryOperator

Symbolic representation of a ? b : c
"""
mutable struct TernaryOperator <: Operator
	name::String
	save::Bool
	a::BooleanExpression
	b::Expression
	c::Expression
	depi::Integer
	depj::Integer
end
#TODO check with true if else cause i can't seem to find a way to override the ternary operator ?
eval_expr(expr::TernaryOperator, vals::AbstractDict) = eval_expr(expr.a, vals) ? eval_expr(expr.b, vals) : eval_expr(expr.c, vals)
string(expr::TernaryOperator) = "vifelse($(string(expr.a[expr.depi, expr.depj])), $(string(expr.b[expr.depi, expr.depj])), $(string(expr.c[expr.depi, expr.depj])))"
prec(expr::TernaryOperator) = 10
getindex(expr::TernaryOperator, depi, depj) = TernaryOperator(expr.name, expr.save, expr.a, expr.b, expr.c, expr.depi+depi, expr.depj+depj)
TernaryOperator(a::BooleanExpression, b::Expression, c::Expression; name="TA_"*a.name*"_"*b.name*"_"*c.name, save=false) = TernaryOperator(name, save, a, b, c, 0, 0)

#Conditionals
"""
	GreaterThan

tests if left > right
"""
mutable struct GreaterThan <: BinaryBooleanOperator
	name::String
	save::Bool
	left::Expression
	right::Expression
	depi::Integer
	depj::Integer
end
symbol(expr::GreaterThan) = ">"
op(expr::GreaterThan) = >
prec(expr::GreaterThan) = 10
GreaterThan(left, right; name="gt_"*left.name*"_"*right.name, save=false) = GreaterThan(name, save, left, right, 0, 0)
>(left::Expression, right::Expression) = GreaterThan(left, right)
>(left::Expression, right::Real) = left > RealValue(right)
>(left::Real, right::Expression) = RealValue(left) > right

