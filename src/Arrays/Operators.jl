"""
Operators
"""
abstract type Operator <: Expression end

"""
Unary Operators
"""
abstract type UnaryOperator <: Operator end
getindex(expr::UnaryOperator, depx, depy) = typeof(expr)(expr.name, expr.expr, expr.depx+depx, expr.depy+depy)
string(expr::UnaryOperator) = "($(symbol(expr))($(string(expr.expr[expr.depx, expr.depy]))))"
eval(expr::UnaryOperator, vals::AbstractDict) = op(expr)(eval(expr.expr, vals))

"""
Negative, represents the negation of an expression
"""
mutable struct Negative <: UnaryOperator
	name::String
	expr::Expression
        depx :: Integer
        depy :: Integer
end
symbol(expr::Negative) = "-"
op(expr::Negative) = -
prec(expr::Negative) = 2
Negative(name, expr) = Negative(name, expr, 0, 0)
-(expr::Expression) = Negative(expr.name*"_neg", expr)

"""
AbsoluteValue
"""
mutable struct AbsoluteValue <: UnaryOperator
	name::String
	expr::Expression
	depx::Integer
	depy::Integer
end
symbol(expr::AbsoluteValue) = "abs"
op(expr::AbsoluteValue) = abs
prec(expr::AbsoluteValue) = 10
AbsoluteValue(name, expr) = AbsoluteValue(name, expr, 0, 0)
abs(expr::Expression) = AbsoluteValue("abs_"*expr.name, expr)

"""
Binary Operators
"""
abstract type BinaryOperator <: Operator end
getindex(expr::BinaryOperator, depx, depy) = typeof(expr)(expr.name, expr.left, expr.right, expr.depx+depx, expr.depy+depy)
string(expr::BinaryOperator) = "($(string(expr.left[expr.depx, expr.depy])))$(symbol(expr))($(string(expr.right[expr.depx, expr.depy])))"
eval(expr::BinaryOperator, vals::AbstractDict) = op(expr)(eval(expr.left, vals), eval(expr.right, vals))

"""
Addition
"""
mutable struct Addition <: BinaryOperator
	name::String
	left::Expression
	right::Expression
        depx :: Integer
        depy :: Integer
end
symbol(expr::Addition) = "+"
op(expr::Addition) = +
prec(expr::Addition) = 1
Addition(name, left, right) = Addition(name, left, right, 0,0)
+(left::Expression, right::Expression) = Addition("p_"*left.name*"_"*right.name, left, right)
+(left::Expression, right::Real) = left + RealValue(right)
+(left::Real, right::Expression) = RealValue(left) + right

"""
Substraction
"""
mutable struct Substraction <: BinaryOperator
	name::String
	left::Expression
	right::Expression
        depx :: Integer
        depy :: Integer
end
symbol(expr::Substraction) = "-"
op(expr::Substraction) = -
prec(expr::Substraction) = 2
Substraction(name, left, right) = Substraction(name, left, right, 0,0)
-(left::Expression, right::Expression) = Substraction("m_"*left.name*"_"*right.name, left, right)
-(left::Expression, right::Real) = Substraction(left, RealValue(right))
-(left::Real, right::Expression) = Substraction(RealValue(left), right)

"""
Multiplication
"""
mutable struct Multiplication <: BinaryOperator
	name::String
	left::Expression
	right::Expression
        depx :: Integer
        depy :: Integer
end
symbol(expr::Multiplication) = "*"
op(expr::Multiplication) = *
prec(expr::Multiplication) = 3
Multiplication(name, left, right) = Multiplication(name, left, right, 0,0)
*(left::Expression, right::Expression) = Multiplication("t_"*left.name*"_"*right.name, left, right)
*(left::Expression, right::Real) = left * RealValue(right)
*(left::Real, right::Expression) = RealValue(left) * right

"""
Division

TODO error handling
"""
mutable struct Division <: BinaryOperator
	name::String
	left::Expression
	right::Expression
        depx :: Integer
        depy :: Integer
end
symbol(expr::Division) = "/"
op(expr::Division) = /
prec(expr::Division) = 3
Division(name, left, right) = Division(name, left, right, 0,0)
/(left::Expression, right::Expression) = Division("d_"*left.name*"_"*right.name, left, right)
/(left::Expression, right::Real) = left / RealValue(right)
/(left::Real, right::Expression) = RealValue(left) / right

"""
Exponentiation
"""
mutable struct Exponentiation <: BinaryOperator
	name::String
	left::Expression
	right::Expression
	depx::Integer
	depy::Integer
end
symbol(expr::Exponentiation) = "^"
op(expr::Exponentiation) = ^
prec(expr::Exponentiation) = 10
Exponentiation(name, left, right) = Exponentiation(name, left, right, 0,0)
^(left::Expression, right::Expression) = Exponentiation("pow_"*left.name*"_"*right.name, left, right)
^(left::Expression, right::Real) = left ^ RealValue(right)
^(left::Real, right::Expression) = RealValue(left) ^ right

abstract type BinaryBooleanOperator <: BinaryOperator end
abstract type UnaryBooleanOperator <: UnaryOperator end

BooleanExpression = Union{UnaryBooleanOperator, BinaryBooleanOperator}

export TernaryOperator
"""
TernaryOperator

	symbolic representation of a ? b : c
"""
mutable struct TernaryOperator <: Operator
	name::String
	a::BooleanExpression
	b::Expression
	c::Expression
	depx::Integer
	depy::Integer
end
#TODO check with true if else cause i can't seem to find a way to override the ternary operator ?
eval(expr::TernaryOperator, vals::AbstractDict) = eval(expr.a, vals) ? eval(expr.b, vals) : eval(expr.c, vals)
string(expr::TernaryOperator) = "($(string(expr.a[expr.depx, expr.depy]))) ? ($(string(expr.b[expr.depx, expr.depy]))) : ($(string(expr.c[expr.depx, expr.depy])))"
prec(expr::TernaryOperator) = 10
getindex(expr::TernaryOperator, depx, depy) = TernaryOperator(expr.name, expr.a, expr.b, expr.c, expr.depx+depx, expr.depy+depy)
TernaryOperator(name, a, b, c) = TernaryOperator(name, a, b, c, 0, 0)
TernaryOperator(a, b, c) = TernaryOperator("TA_" * a.name * "_" * b.name * "_" * c.name, a, b, c)

#Conditionals
"""
GreaterThan

	tests if left > right
"""
mutable struct GreaterThan <: BinaryBooleanOperator
	name::String
	left::Expression
	right::Expression
	depx::Integer
	depy::Integer
end
symbol(expr::GreaterThan) = ">"
op(expr::GreaterThan) = >
prec(expr::GreaterThan) = 10
GreaterThan(name, left, right) = GreaterThan(name, left, right, 0, 0)
>(left::Expression, right::Expression) = GreaterThan(left.name*"_"*right.name*"_gt", left, right)
>(left::Expression, right::Real) = left > RealValue(right)
>(left::Real, right::Expression) = RealValue(left) > right

