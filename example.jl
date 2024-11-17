using Revise
using glOpqua
using Plots
using LinearAlgebra
using BenchmarkTools

n_strains = 7

#TODO: add antigenic mapping that based on distance map (or pairwise distances),

# computes pairwise distances (if necessary), cross-immunities, and imprinting and maturation fractions.
immunity_coef = [
    # This tracks strain cross immunity (self immunity is on diagonal).
    # Row number (first index) is pre-existing immunity to that number strain,
    # columns are number of infecting strain.
    1.0 0.0
    0.0 1.0
]
# or
immunity_coef = Matrix{Float64}(I, n_strains, n_strains)

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
parameters = glOpqua.parameters(
    n_strains=n_strains,
    transmission_rates=repeat([0.1], outer=n_strains),
    recovery_rates=repeat([0.05], outer=n_strains),
    imprinted_loss_rate=0.01,
    matured_loss_rate=0.1,
    frac_imprinted=frac_imprinted,
    frac_matured=frac_matured,
    immunity_coef=immunity_coef,
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
end

# Time
tspan = (0.0, 200.0)

# Simulation
ode_sol = glOpqua.simulate(parameters, init, tspan)

@benchmark ode_sol_bm = glOpqua.simulate(parameters, init, tspan)

# Data processing
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
