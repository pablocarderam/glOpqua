n_strains = 2

n_compartments = n_strains + 2#2 * 3^n_strains

transmission_rate = 0.1
recovery_rate = 0.1

using DifferentialEquations
using Plots

function sir_ode2(du,u,p,t)
    #S,I,R = u
    #b,g = p
    du[1] = -p[1]*u[1]*u[2]
    du[2] = p[1]*u[1]*u[2]-p[2]*u[2]
    du[3] = p[2]*u[2]
end
parms = [0.1,0.05]
init = [0.99,0.01,0.0]
tspan = (0.0,200.0)
sir_prob2 = ODEProblem(sir_ode2,init,tspan,parms)
sir_sol = solve(sir_prob2,saveat = 0.1)

plot(sir_sol,xlabel="Time",ylabel="Number")

@benchmark sir_sol = solve(sir_prob2,QNDF(),saveat = 0.1)
