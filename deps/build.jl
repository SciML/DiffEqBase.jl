infomsg = """
~~~ DifferentialEquations.jl *BREAKING* changes  ~~~

We have changed the front-end API on how
users may define equations of motion and 
problems, for all problem types that can
be used in the DifferentialEquations.jl.
These are *BREAKING* changes, and they
also have *NO WARNINGS*!

Please see our latest documentation here:
http://docs.juliadiffeq.org/latest/

or the blogpost that describes the changes:
http://juliadiffeq.org/2018/01/24/Parameters.html

In short, the mutated argument is the first argument,
and parameters are now directly passed
into the equations of motion function. For all
types now mutation goes first, then dependent variables, 
then parameters, then independent variables. 

`f(mutated, dependent variables, p/integrator, independent variables)`

For example, this means that the ODE syntax will be `f(u,p,t)` (for the
out-of-place) and `f(du,u,p,t)` (for the in-place). Notice
that this change also removes the need for ParameterizedFunctions
as now parameters are part of the equations of motion.

For more details please visit the above links!
"""

info(infomsg)
