using StaticArrays
using Combinatorics
using DifferentialEquations
using Plots

# const NAIVE::Int64 = 1
# const IMPRINTED::Int64 = 2
# const MATURED::Int64 = 3

const NAIVE::Int64 = 0
const IMPRINTED::Int64 = 1
const MATURED::Int64 = 2

const NAIVE_VEC::Vector{Bool} = [1, 0, 0]
const IMPRINTED_VEC::Vector{Bool} = [0, 1, 0]
const MATURED_VEC::Vector{Bool} = [0, 0, 1]

n_strains = 2
n_immunities = 3^n_strains # number of types of compartments, i.e. possible immune state combinations
n_compartments = (n_strains + 1) * n_immunities # n_strains possible infecting strains plus uninfected per each immune compartment type

@kwdef struct ModelParameters
    n_strains::Int64
    transmission_rate::Float64
    recovery_rate::Float64
end

params = ModelParameters(
    transmission_rate=0.1,
    recovery_rate=0.05,
    frac_imprinted=1.0, # this will be a map
    frac_matured=0.0, # this will be a map
    n_strains=2
)

immunities_ids = reduce(vcat, collect(multiset_permutations(i, n_strains)) for i in with_replacement_combinations([NAIVE, IMPRINTED, MATURED], n_strains))
immunities_dict = Dict{Vector{Int64},Int64}(zip(immunities_ids, collect(1:n_immunities)))
immunities_list = reduce(vcat, collect(multiset_permutations(i, n_strains)) for i in with_replacement_combinations([NAIVE_VEC, IMPRINTED_VEC, MATURED_VEC], n_strains))
immunities = cat(dims=3, [hcat(compartment...) for compartment in immunities_list]...)
# Dimensions of this matrix are: immune state (naive, imprinted, matured):strains:host compartment
# This entire matrix appliers to both uninfected and infected with each strain, so there are (1 + n_strains) as many compartments as dimension 3 of this matrix

function odeFunc(du, u, p::ModelParameters, t::Float64)
    # du[1] = -p.transmission_rate * u[1] * sum(@views u[2:p.n_strains+1])

    # Temp variables for the loop:
    infected = zeros(n_immunities)
    naive = zeros(n_strains)
    imprinted = zeros(n_strains)
    matured = zeros(n_strains)
    immunity = zeros(n_strains)
    for c in 1:n_immunities
        naive = zeros(n_strains)
        imprinted = zeros(n_strains)
        matured = zeros(n_strains)
        # Calculate totals:
        for s in 1:n_strains
            infected[s] = p.transmission_rate * u[c] * sum(@views u[(s*n_immunities+1):((s+1)*n_immunities)])

            immunity = immunities_ids[c]
            immunity[s] = 0 # We consider the immune status identical to the current one but with naive in the current strain
            # I think immunities_dict[immunities_ids[c]] == c, so instead of u[immunities_dict[immunities_ids[c]]+(s*n_immunities)] in these next three blocks, we do u[c+(s*n_immunities)]
            if immunities_ids[c][s] == IMPRINTED
                imprinted[s] = p.recovery_rate * (p.frac_imprinted * u[immunities_dict[immunity]+(s*n_immunities)] + u[c+(s*n_immunities)])
            elseif immunities_ids[c][s] == MATURED
                matured[s] = p.recovery_rate * (p.frac_matured * u[immunities_dict[immunity]+(s*n_immunities)] + u[c+(s*n_immunities)])
            else
                naive[s] = p.recovery_rate * ((1 - p.frac_imprinted - p.frac_matured) * u[immunities_dict[immunity]+(s*n_immunities)] + u[c+(s*n_immunities)])
            end
        end

        # Uninfected compartment for this immunity type:
        du[c] = -sum(infected) + sum(imprinted) + sum(matured) + sum(naive)

        # Infected compartments for this immunity type:
        for s in 1:n_strains
            du[c+(s*n_immunities)] = infected[s] - p.recovery_rate * 1 * u[c+(s*n_immunities)]
        end
    end

    # # Infected compartments
    # for s in 1:n_strains
    #     for c in 1:n_immunities
    #         # Compartment order:
    #         # 1 through 2^(n_strains): Uninfected with immunity to
    #         du[i] = p.transmission_rate * sum(u[1]) * u[i] - p.recovery_rate * u[i]
    #         du[p.n_strains+i] = p.recovery_rate * u[i] #p.transmission_rate * u[1] * u[i] - p.recovery_rate * u[i]
    #     end
    # end
    # # du[p.n_strains+2] = p.recovery_rate * sum(@views u[2:p.n_strains+1])
end

# transmission_rate = 0.1
# recovery_rate = 0.05

# params = [transmission_rate, recovery_rate, n_strains]

start_freq = 0.01
init = zeros(n_compartments)
init[1] = 1.0 - (start_freq * n_strains)
for i in 1:n_strains
    init[i*n_immunities+1] = start_freq
end

tspan = (0.0, 200.0)
ode_prob = ODEProblem(odeFunc, init, tspan, params)
ode_sol = solve(ode_prob, QNDF(), saveat=0.1);

plot(ode_sol, xlabel="Time", ylabel="Number", linewidth=2)

# @benchmark ode_sol = solve(odeFunc, QNDF(), saveat=0.1)
