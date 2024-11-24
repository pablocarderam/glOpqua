using Revise
using glOpqua
using Plots
using LinearAlgebra
using BenchmarkTools

n_strains = 6
# 7 gives a stackoverflow error if using default solver option, but runs on QNDF(autodiff=false).
# For n_strains < 7, the default solver is faster than QNDF.

# Computes pairwise distances (if necessary), evolutionary risk motes, cross-immunities, and imprinting and maturation fractions.
#TODO: add antigenic mapping that based on distance map (or pairwise distances),

immunity_coef = [
    # This tracks strain cross immunity (self immunity is on diagonal).
    # Row number (first index) is pre-existing immunity to that number strain,
    # columns are number of infecting strain.
    1.0 0.0
    0.0 1.0
]
# or
immunity_coef = Matrix{Float64}(I, n_strains, n_strains)

evorisk_coef = ones(2, 2, 3)
evorisk_coef[:, :, 1] = [
    # This tracks evolutionary risk motes for naive immunity.
    # Should be a diagonal matrix, one value per strain.
    0.01 0.0
    0.0 0.01
]
evorisk_coef[:, :, glOpqua.IMPRINTED+1] = [
    # This tracks evolutionary risk motes for imprinted immunity.
    # Row number (first index) is pre-existing immunity to that number strain,
    # columns are number of infecting strain.
    0.02 0.1
    0.1 0.02
]
evorisk_coef[:, :, glOpqua.MATURED+1] = [
    # This tracks evolutionary risk motes for affinity matured immunity.
    # Row number (first index) is pre-existing immunity to that number strain,
    # columns are number of infecting strain.
    0.04 0.2
    0.2 0.04
]
# or
evorisk_coef = ones(n_strains, n_strains, 3)
evorisk_coef[:, :, 1] = 0.0 * Matrix{Float64}(I, n_strains, n_strains)
evorisk_coef[:, :, glOpqua.IMPRINTED+1] = evorisk_coef[:, :, glOpqua.IMPRINTED] + 0.0 * Matrix{Float64}(I, n_strains, n_strains)
evorisk_coef[:, :, glOpqua.MATURED+1] = evorisk_coef[:, :, glOpqua.MATURED] + 0.0 * Matrix{Float64}(I, n_strains, n_strains)

frac_imprinted = [
    # This tracks likelihood of developing an imprinted response to a new strain given
    # pre-existing immunity to another strain.
    # Diagonal is zero (immunity is not re-imprinted on infection with same strain).
    # Row number (first index) is pre-existing immunity to that number strain,
    # columns are number of infecting strain.
    0.0 0.1
    0.1 0.0
]
# or
frac_imprinted = 0.1 * (ones(n_strains, n_strains) - Matrix{Float64}(I, n_strains, n_strains))

frac_matured = [
    # This tracks likelihood of developing an affinity matured response to a new strain given
    # pre-existing immunity to another strain.
    # Diagonal is zero (it is true immunity is affinity-matured on reexposure,
    # but negative effects on antigenic escape likelihood are avoided).
    # Row number (first index) is pre-existing immunity to that number strain,
    # columns are number of infecting strain.
    0.0 0.9
    0.9 0.0
]
# or
frac_matured = ones(n_strains, n_strains) - frac_imprinted - Matrix{Float64}(I, n_strains, n_strains)
# Notice frac_imprinted and frac_matured must sum to <=1 at each element.
# 1-frac_imprinted-frac_matured is the fraction of cases with preexisting immunity
# to the row number strain that do not acquire any immunity to the column number
# strain after recovering from that infection (this diagonal is not taken into
# account inside the ODE function).

# Parameters
#TODO: add non-immune compartment-affecting strains
#TODO: add option for single immune compartment type and just using evorisk to track impact of imprinted vs matured??? nah this doesn't work
parameters = glOpqua.parameters(
    n_strains=n_strains,
    transmission_rates=repeat([0.1], outer=n_strains),
    recovery_rates=repeat([0.05], outer=n_strains),
    imprinted_loss_rate=0.01,
    matured_loss_rate=0.1,
    frac_imprinted=frac_imprinted,
    frac_matured=frac_matured,
    immunity_coef=immunity_coef,
    evorisk_coef=evorisk_coef,
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

plot(
    ode_sol,
    xlabel="Time", ylabel="Number", linewidth=2,
    # labels=reshape(vcat(["Uninfected"],["Strain "*string(i) for i in 1:parameters.n_strains]), 1, 1+parameters.n_strains),
)
