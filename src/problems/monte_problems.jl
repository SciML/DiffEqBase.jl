type MonteCarloProblem{T,T2,T3} <: AbstractMonteCarloProblem
  prob::T
  prob_func::T2
  output_func::T3
end

MonteCarloProblem(prob::DEProblem;
                  output_func = (sol,i)-> sol,
                  prob_func= (prob,i)->prob) =
                  MonteCarloProblem(prob,prob_func,output_func)
