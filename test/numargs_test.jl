using DiffEqBase

function test_num_args()
  f(x) = 2x
  f(x,y) = 2xy

  numpar = DiffEqBase.numargs(f) # Should be 2
  g = (x,y) -> x^2
  numpar2 = DiffEqBase.numargs(g)
  numpar == 2 && numpar2 == 2
end
@test test_num_args()
