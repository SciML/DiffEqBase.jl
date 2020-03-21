using DiffEqBase: bisection
using Test

# basic tests
# result of bisection should never be zero, 
# should have sign same as left limit
# and nexfloat in tdir should be zero or same sign of right limit
r = bisection(sin, (3.0, 4.0), 1.0)
@test !iszero(sin(r))
@test sin(r) > 0
@test sin(nextfloat(r)) <= 0

r = bisection(sin, (4.0, 3.0), -1.0)
@test !iszero(sin(r))
@test sin(r) < 0
@test sin(prevfloat(r)) >= 0

# this example has a big zero interval, makes it go inside the inner loop
for i in 1:10
	f = (x) -> (10.0^(-i))*(x^3)
	r = bisection(f, (-1.0, 1.0), 1.0)
	@test !iszero(f(r))
	@test f(r) < 0
	@test f(nextfloat(r)) >= 0

	r = bisection(f, (1.0, -1.0), -1.0)
	@test !iszero(f(r))
	@test f(r) > 0
	@test f(prevfloat(r)) <= 0
end
