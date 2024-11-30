using Revise
using glOpqua
using Plots
using LinearAlgebra
using BenchmarkTools

n_ag_space_dimensions = 2
n_gen_space_dimensions = 2
# The number of dimensions of antigenic and genetic space does not need to be equal.

n_ag_clusters = 5
# =8 takes 13.78 s, =9 takes 130.99 s to solve for time=200

n_strains = ones(Int, n_ag_clusters)#[1,2,2,1]
# How many strains are in each antigenic cluster

ag_cluster_positions = zeros(n_ag_clusters, n_ag_space_dimensions)
ag_cluster_positions[:, end] .+= 10.0 .* collect(1:n_ag_clusters)

strain_positions = zeros(sum(n_strains), n_gen_space_dimensions)
strain_positions[:, end] .+= 10.0 .* collect(1:sum(n_strains))

# Compute pairwise distances, cross-immunities, and imprinting and maturation fractions.
ag_map = glOpqua.antigenicMap(
    positions=ag_cluster_positions,
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
    K_evorisk=0.5 .* ones(glOpqua.N_IMMUNE_STATES),
    h_evorisk=5.0 .* ones(glOpqua.N_IMMUNE_STATES),
    naive_evorisk=ones(sum(n_strains)) .* 0.0,
    imprinted_matured_evorisk=ones(sum(n_strains), glOpqua.N_IMMUNE_STATES - 1) .* 1.0,)

# Parameters
#TODO: add option for single immune compartment type and just using evorisk to track impact of imprinted vs matured??? nah this doesn't work
parameters = glOpqua.parameters(
    n_ag_clusters=n_ag_clusters,
    n_strains=n_strains,
    transmission_rates=[repeat([0.1], outer=n_strains[i]) for i in 1:n_ag_clusters],
    recovery_rates=[repeat([0.05], outer=n_strains[i]) for i in 1:n_ag_clusters],
    imprinted_loss_rate=0.01,
    matured_loss_rate=0.1,
    frac_imprinted=ag_map.frac_imprinted,
    frac_matured=ag_map.frac_matured,
    immunity_coef=ag_map.immunity_coef,
    evorisk_coef=evo_map.evorisk_coef,
    birth_rate_uninfected=0.01,
    birth_rates_infected=[repeat([0.01], outer=n_strains[i]) for i in 1:n_ag_clusters],
    death_rate_uninfected=0.01,
    death_rates_infected=[repeat([0.01], outer=n_strains[i]) for i in 1:n_ag_clusters],)

# Initial conditions
start_freq = 0.01
init = zeros(parameters.n_compartments)
init[1] = 1.0 - (start_freq * sum(parameters.n_strains))
for i in 1:sum(parameters.n_strains)
    init[i*parameters.n_immunities+1] = start_freq
    init[(sum(parameters.n_strains)+1)*parameters.n_immunities+i] = start_freq * parameters.evorisk_coef[(repeat([glOpqua.NAIVE], outer=parameters.n_ag_clusters), i)]
end

# Time
tspan = (0.0, 200.0)

# Simulation
#TODO: add piecewise ODE with stochastic strain appearance and extinction
ode_sol = glOpqua.simulate(parameters, init, tspan)
# ode_sol = glOpqua.simulate(parameters, init, tspan, solver="QNDF")
# @benchmark ode_sol_bm = glOpqua.simulate(parameters, init, tspan)

# Data processing
dat_evo = glOpqua.evoriskData(ode_sol, parameters)
plot(
    dat_evo[:, :t],
    vcat([dat_evo[:, "Strain_"*string(i)] for i in 1:sum(parameters.n_strains)]),
    xlabel="Time", ylabel="Number", linewidth=2,
    labels=reshape(vcat(["Strain " * string(i) for i in 1:sum(parameters.n_strains)]), 1, sum(parameters.n_strains)),
)

dat_agg = glOpqua.dataByStrain(ode_sol, parameters)
plot(
    dat_agg[:, :t],
    vcat([dat_agg[:, "Uninfected"]], [dat_agg[:, "Strain_"*string(i)] for i in 1:sum(parameters.n_strains)]),
    xlabel="Time", ylabel="Number", linewidth=2,
    labels=reshape(vcat(["Uninfected"], ["Strain " * string(i) for i in 1:sum(parameters.n_strains)]), 1, 1 + sum(parameters.n_strains)),
)

# plot(
#     ode_sol,
#     xlabel="Time", ylabel="Number", linewidth=2,
#     # labels=reshape(vcat(["Uninfected"],["Strain "*string(i) for i in 1:parameters.n_ag_clusters]), 1, 1+parameters.n_ag_clusters),
# )
