abstract type UnaryBooleanOperator <: UnaryOperator end

abstract type BinaryBooleanOperator <: BinaryOperator end

BooleanExpression = Union{UnaryBooleanOperator, BinaryBooleanOperator}

"""
GreaterThan

	tests if left > right
"""
struct GreaterThan{L<:Expression, R<:Expression} <: BinaryBooleanOperator
	left::L
	right::R
end
symbol(expr::GreaterThan) = ">"
op(expr::GreaterThan) = >
prec(expr::GreaterThan) = 10
>(left::Expression, right::Expression) = GreaterThan(left, right)
>(left::Expression, right::Real) = GreaterThan(left, RealValue(right))
>(left::Real, right::Expression) = GreaterThan(RealValue(left), right)


export TernaryOperator
"""
TernaryOperator

	symbolic representation of a ? b : c
"""
struct TernaryOperator{L<:Expression, R<:Expression} <: Operator
	a::BooleanExpression
	b::L
	c::R
end
#TODO check with true if else cause i can't seem to find a way to override the ternary operator ?
eval(expr::TernaryOperator, vals::AbstractDict) = eval(expr.a, vals) ? eval(expr.b, vals) : eval(expr.c, vals)
string(expr::TernaryOperator) = "($(string(expr.a))) ? ($(string(expr.b))) : ($(string(expr.c)))"
	
