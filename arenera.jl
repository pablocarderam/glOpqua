using StaticArrays
using DifferentialEquations
using Plots

n_strains = 2

n_compartments = n_strains + 2#2 * 3^n_strains

function odeFunc(du,u,p,t)
    #S,I,R = u
    #b,g = p
    du[1] = -p[1]*u[1]*sum(@views u[2:p[3]+1])
    for i in 2:p[3]+1
        du[i] = p[1]*u[1]*u[i]-p[2]*u[i]
    end
    du[p[3]+2] = p[2]*sum(@views u[2:p[3]+1])
end

transmission_rate = 0.1
recovery_rate = 0.05

params = SA[transmission_rate, recovery_rate, n_strains]

init = [0.99,0.01,0.0]
tspan = (0.0,200.0)
ode_prob = ODEProblem(odeFunc,init,tspan,params)
ode_sol = solve(ode_prob,QNDF(),saveat = 0.1)

plot(ode_sol,xlabel="Time",ylabel="Number",linewidth=2)

@benchmark ode_sol = solve(odeFunc,QNDF(),saveat = 0.1)
