using Revise
using glOpqua
using Plots
using BenchmarkTools

n_strains = 2

immunity_coef = [
    # rows are immunities, columns are infecting strain
    [1.0,0.0],
    [0.0,1.0]
    ]

# Parameters
params = glOpqua.parameters(
    n_strains=n_strains,
    transmission_rates=repeat([0.1],outer=n_strains),
    recovery_rates=repeat([0.05],outer=n_strains),
    frac_imprinted=0.5, # this will be a map
    frac_matured=0.5, # this will be a map
    immunity_coef=immunity_coef,
)

# Initial conditions
start_freq = 0.01
init = zeros(params.n_compartments)
init[1] = 1.0 - (start_freq * params.n_strains)
for i in 1:params.n_strains
    init[i*params.n_immunities+1] = start_freq
end

# Time
tspan = (0.0, 200.0)

ode_sol = glOpqua.simulate(params,init,tspan)
dat_agg = glOpqua.dataByStrain(ode_sol,params)
plot(
    dat_agg[:, :t],
    vcat([dat_agg[:, "Uninfected"]], [dat_agg[:, "Strain_"*string(i)] for i in 1:params.n_strains]),
    xlabel="Time", ylabel="Number", linewidth=2
)
