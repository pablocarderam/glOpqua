using StaticArrays
using DifferentialEquations
using Plots

n_strains = 20

n_compartments = n_strains + 2#2 * 3^n_strains

# @kwdef struct ModelParameters
#     n_strains::Int64
#     transmission_rate::Float64
#     recovery_rate::Float64
# end

# params = ModelParameters(
#     transmission_rate=0.1,
#     recovery_rate=0.05,
#     n_strains=2
# )

function odeFunc(du, u, p, t)
    #S,I,R = u
    #b,g = p
    du[1] = -p.transmission_rate * u[1] * sum(@views u[2:p.n_strains+1])
    for i in 2:p.n_strains+1
        du[i] = p.transmission_rate * u[1] * u[i] - p.recovery_rate * u[i]
        # du[i] = p.transmission_rate * u[1] * u[i] - p.recovery_rate * u[i]
    end
    du[p.n_strains+2] = p.recovery_rate * sum(@views u[2:p.n_strains+1])
end

transmission_rate = 0.1
recovery_rate = 0.05

params = SA[transmission_rate, recovery_rate, n_strains]

function odeFunc2(du::Vector{Float64}, u::Vector{Float64}, p::Vector{Float64}, t::Float64)
    #S,I,R = u
    #b,g = p
    du[1] = -p[1] * u[1] * sum(@views u[2:floor(Int, p[3])+1])
    for i in 2:floor(Int, p[3])+1
        du[i] = p[1] * u[1] * u[i] - p[2] * u[i]
        du[i] = p[1] * u[1] * u[i] - p[2] * u[i]
    end
    du[floor(Int, p[3])+2] = p[2] * sum(@views u[2:floor(Int, p[3])+1])
end

start_freq = 0.01
init = vcat([1.0 - (start_freq * n_strains)], start_freq * ones(n_strains), [0.0])
tspan = (0.0, 200.0)
ode_prob = ODEProblem(odeFunc2, init, tspan, params)
ode_sol = solve(ode_prob, QNDF(), saveat=0.1);

plot(ode_sol, xlabel="Time", ylabel="Number", linewidth=2)

@benchmark ode_sol = solve(odeFunc, QNDF(), saveat=0.1)
