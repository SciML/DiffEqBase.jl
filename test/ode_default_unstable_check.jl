using Test, RecursiveArrayTools, StaticArrays

using DiffEqBase: NAN_CHECK

@test NAN_CHECK(3.0+4.0im) == false
@test NAN_CHECK(NaN) == true

u1 = ones(3)
@test NAN_CHECK(u1) == false
u1′ = copy(u1)
u1′[2] = NaN
@test NAN_CHECK(u1′) == true


u2 = [SA[1.0 1.0; 1.0 1.0] for i = 1:3]
@test NAN_CHECK(u2) == false
u2′ = copy(u2)
u2′[2] = SA[1.0 NaN; 1.0 1.0]
@test NAN_CHECK(u2′) == true

u3 = VectorOfArray([ones(5), ones(5)])
@test NAN_CHECK(u3) == false
u3′ = recursivecopy(u3)
u3′[2][3] = NaN
@test NAN_CHECK(u3′) == true

u4 = ArrayPartition(u1, u2, u3)
@test NAN_CHECK(u4) == false
u4_1 = ArrayPartition(u1′, u2, u3)
@test NAN_CHECK(u4_1) == true
u4_2 = ArrayPartition(u1, u2′, u3)
@test NAN_CHECK(u4_2) == true
u4_3 = ArrayPartition(u1, u2, u3′)
@test NAN_CHECK(u4_3) == true

@test NAN_CHECK(ArrayPartition(u4, u4)) == false
@test NAN_CHECK(ArrayPartition(u4, u4_1)) == true
@test NAN_CHECK(ArrayPartition(u4, u4_2)) == true
@test NAN_CHECK(ArrayPartition(u4, u4_3)) == true
