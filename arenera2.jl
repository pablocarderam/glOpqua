using Revise
using BenchmarkTools
using StaticArrays
using Combinatorics
using DifferentialEquations
using DataFrames
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

@kwdef struct ModelParameters
    n_strains::Int64
    transmission_rate::Float64
    recovery_rate::Float64
    frac_imprinted::Float64
    frac_matured::Float64
end

@kwdef struct ModelParametersInner
    n_strains::Int64
    n_immunities::Int64
    n_compartments::Int64
    immunities_ids::Vector{Vector{Int64}}
    immunities_dict::Dict{Vector{Int64},Int64}
    immunities_list::Vector{Vector{Vector{Bool}}}
    immunities::Array{Bool, 3} # 3 comes from 3 immune states, hardcoded
    transmission_rate::Float64
    recovery_rate::Float64
    frac_imprinted::Float64
    frac_matured::Float64
end

function simulate(params::ModelParameters, odeFunc::Function)
    n_strains = params.n_strains
    n_immunities = 3^n_strains # number of types of compartments, i.e. possible immune state combinations
    n_compartments = (n_strains + 1) * n_immunities # n_strains possible infecting strains plus uninfected per each immune compartment type

    immunities_ids = reduce(vcat, collect(multiset_permutations(i, n_strains)) for i in with_replacement_combinations([NAIVE, IMPRINTED, MATURED], n_strains))
    immunities_dict = Dict{Vector{Int64},Int64}(zip(immunities_ids, collect(1:n_immunities)))
    immunities_list = reduce(vcat, collect(multiset_permutations(i, n_strains)) for i in with_replacement_combinations([NAIVE_VEC, IMPRINTED_VEC, MATURED_VEC], n_strains))
    immunities = cat(dims=3, [hcat(compartment...) for compartment in immunities_list]...)
    # Dimensions of this matrix are: immune state (naive, imprinted, matured):strains:host compartment
    # This entire matrix appliers to both uninfected and infected with each strain, so there are (1 + n_strains) as many compartments as dimension 3 of this matrix
    paramsinner = ModelParametersInner(
        n_strains=n_strains,
        n_immunities=n_immunities,
        n_compartments=n_compartments,
        immunities_ids=immunities_ids,
        immunities_dict=immunities_dict,
        immunities_list=immunities_list,
        immunities=immunities,
        transmission_rate=params.transmission_rate,
        recovery_rate=params.recovery_rate,
        frac_imprinted=params.frac_imprinted,
        frac_matured=params.frac_matured,
    )

    start_freq = 0.01
    init = zeros(n_compartments)
    init[1] = 1.0 - (start_freq * n_strains)
    for i in 1:n_strains
        init[i*n_immunities+1] = start_freq
    end

    tspan = (0.0, 200.0)
    ode_prob = ODEProblem(odeFunc, init, tspan, paramsinner)
    ode_sol = solve(ode_prob, QNDF(autodiff=false), saveat=0.1)

    dat = DataFrame(Tables.table(ode_sol'))
    dat[!, :t] = ode_sol.t

    dat_agg = DataFrame()
    dat_agg[!, :t] = ode_sol.t

    dat_agg[!, ("Uninfected")] = sum(eachcol(dat[:, 1:n_immunities]))
    for i in 1:n_strains
        dat_agg[!, ("Strain_"*string(i))] = sum(eachcol(dat[:, (n_immunities*i+1):(n_immunities*(i+1))]))
    end

    return dat_agg
end

function odeFunc(du, u, p::ModelParametersInner, t::Float64)
    # Temp variables for the loop:
    infected::Vector{Float64} = zeros(Float64,p.n_immunities)
    immunity::Vector{Int64} = zeros(Int64,p.n_strains)

    for c in 1:p.n_immunities
        naive = 0.0
        imprinted = 0.0
        matured = 0.0
        total_infected = 0
        # Calculate totals:
        for s in 1:p.n_strains
            infected[s] = p.transmission_rate * u[c] * sum(@views u[(s*p.n_immunities+1):((s+1)*p.n_immunities)])
            total_infected += infected[s]
            # It should be true that immunities_dict[immunities_ids[c]] == c, so instead of u[immunities_dict[immunities_ids[c]]+(s*n_immunities)] in these next three blocks, we do u[c+(s*n_immunities)]
            if p.immunities_ids[c][s] == IMPRINTED
                immunity = copy(p.immunities_ids[c])
                immunity[s] = 0 # We consider the immune status identical to the current one but with naive in the current strain
                imprinted += p.recovery_rate * (p.frac_imprinted * u[p.immunities_dict[immunity]+(s*p.n_immunities)] + u[c+(s*p.n_immunities)])
            elseif p.immunities_ids[c][s] == MATURED
                immunity = copy(p.immunities_ids[c])
                immunity[s] = 0 # We consider the immune status identical to the current one but with naive in the current strain
                matured += p.recovery_rate * (p.frac_matured * u[p.immunities_dict[immunity]+(s*p.n_immunities)] + u[c+(s*p.n_immunities)])
            elseif p.immunities_ids[c][s] == NAIVE
                naive += p.recovery_rate * (1 - p.frac_imprinted - p.frac_matured) * u[c+(s*p.n_immunities)]
            end
        end

        # Uninfected compartment for this immunity type:
        du[c] = -total_infected + imprinted + matured + naive

        # Infected compartments for this immunity type:
        for s in 1:p.n_strains
            du[c+(s*p.n_immunities)] = infected[s] - p.recovery_rate * u[c+(s*p.n_immunities)]
        end
    end
end

params = ModelParameters(
    n_strains=6,
    transmission_rate=0.1,
    recovery_rate=0.05,
    frac_imprinted=1.0, # this will be a map
    frac_matured=0.0, # this will be a map
)

dat_agg = simulate(params,odeFunc)
plot(dat_agg[:, :t], vcat([dat_agg[:, "Uninfected"]], [dat_agg[:, "Strain_"*string(i)] for i in 1:params.n_strains]), xlabel="Time", ylabel="Number", linewidth=2)
# plot(ode_sol, xlabel="Time", ylabel="Number", linewidth=2)
@benchmark ode_sol = solve(odeFunc, QNDF(), saveat=0.1)
