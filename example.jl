using Revise
using glOpqua
using Plots
using LinearAlgebra
using BenchmarkTools

n_space_dimensions = 2
# In this case, used for antigenic and genetic space, but these can have different dimensions.

n_strains = 3
# 7 gives a stackoverflow error if using default solver option, but runs on QNDF(autodiff=false).
# For n_strains < 7, the default solver is faster than QNDF.

strain_positions = zeros(n_strains, n_space_dimensions)
strain_positions[:, end] .+= 10.0 .* collect(1:n_strains)

# Compute pairwise distances, cross-immunities, and imprinting and maturation fractions.
ag_map = glOpqua.antigenicMap(
    positions=strain_positions,
    max_cross_immunity=1.0,
    K_cross_immunity=1.0,
    h_cross_immunity=2.0,
    max_imprinted=1.0,
    K_imprinted=5.0,
    h_imprinted=2.0,
    max_matured=1.0,
    K_matured=5.0,
    h_matured=2.0,)

# Evolutionary risk motes
evo_map = glOpqua.evoRiskMap(
    positions=strain_positions,
    K_evorisk=0.5 .* ones(n_strains),
    h_evorisk=5.0 .* ones(n_strains),
    naive_evorisk=ones(n_strains),
    imprinted_matured_evorisk=ones(n_strains, 2),)

# Parameters
#TODO: add non-immune compartment-affecting strains
#TODO: add option for single immune compartment type and just using evorisk to track impact of imprinted vs matured??? nah this doesn't work
parameters = glOpqua.parameters(
    n_strains=n_strains,
    transmission_rates=repeat([0.1], outer=n_strains),
    recovery_rates=repeat([0.05], outer=n_strains),
    imprinted_loss_rate=0.01,
    matured_loss_rate=0.1,
    frac_imprinted=ag_map.frac_imprinted,
    frac_matured=ag_map.frac_matured,
    immunity_coef=ag_map.immunity_coef,
    evorisk_coef=evo_map.evorisk_coef,
    birth_rate_uninfected=0.01,
    birth_rates_infected=repeat([0.01], outer=n_strains),
    death_rate_uninfected=0.01,
    death_rates_infected=repeat([0.01], outer=n_strains),)

# Initial conditions
start_freq = 0.01
init = zeros(parameters.n_compartments)
init[1] = 1.0 - (start_freq * parameters.n_strains)
for i in 1:parameters.n_strains
    init[i*parameters.n_immunities+1] = start_freq
    init[(parameters.n_strains+1)*parameters.n_immunities+i] = start_freq * parameters.evorisk_coef[(repeat([glOpqua.NAIVE], outer=parameters.n_strains), i)]
end

# Time
tspan = (0.0, 2000.0)

# Simulation
ode_sol = glOpqua.simulate(parameters, init, tspan)
#TODO: add piecewise ODE with stochastic strain appearance and extinction

# @benchmark ode_sol_bm = glOpqua.simulate(parameters, init, tspan)

# Data processing
dat_evo = glOpqua.evoriskData(ode_sol, parameters)
plot(
    dat_evo[:, :t],
    vcat([dat_evo[:, "Strain_"*string(i)] for i in 1:parameters.n_strains]),
    xlabel="Time", ylabel="Number", linewidth=2,
    labels=reshape(vcat(["Strain " * string(i) for i in 1:parameters.n_strains]), 1, parameters.n_strains),
)

dat_agg = glOpqua.dataByStrain(ode_sol, parameters)
plot(
    dat_agg[:, :t],
    vcat([dat_agg[:, "Uninfected"]], [dat_agg[:, "Strain_"*string(i)] for i in 1:parameters.n_strains]),
    xlabel="Time", ylabel="Number", linewidth=2,
    labels=reshape(vcat(["Uninfected"], ["Strain " * string(i) for i in 1:parameters.n_strains]), 1, 1 + parameters.n_strains),
)

# plot(
#     ode_sol,
#     xlabel="Time", ylabel="Number", linewidth=2,
#     # labels=reshape(vcat(["Uninfected"],["Strain "*string(i) for i in 1:parameters.n_strains]), 1, 1+parameters.n_strains),
# )
