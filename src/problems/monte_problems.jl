mutable struct MonteCarloProblem{T,T2,T3,T4,T5} <: AbstractMonteCarloProblem
  prob::T
  prob_func::T2
  output_func::T3
  reduction::T4
  u_init::T5
end

MonteCarloProblem(prob::DEProblem;
    output_func = (sol,i)-> (sol,false),
    prob_func= (prob,i,repeat)->prob,
    reduction = (u,data,I)->(append!(u,data),false),
    u_init = []) =
    MonteCarloProblem(prob,prob_func,output_func,reduction,u_init)
