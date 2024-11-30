using Combinatorics

const NAIVE::Int64 = 0
const IMPRINTED::Int64 = 1
const MATURED::Int64 = 2

const NAIVE_VEC::Vector{Bool} = [1, 0, 0]
const IMPRINTED_VEC::Vector{Bool} = [0, 1, 0]
const MATURED_VEC::Vector{Bool} = [0, 0, 1]

function hillFunction(distance::Float64, K::Float64, h::Float64)
    return distance^h / (K^h + distance^h)
end

@kwdef struct ModelParameters
    n_strains::Int64
    n_immunities::Int64
    n_compartments::Int64
    immunities_ids::Vector{Vector{Int64}}
    immunities_dict::Dict{Vector{Int64},Int64}
    immunities_list::Vector{Vector{Vector{Bool}}}
    immunities::Array{Bool,3} # 3 comes from 3 immune states, hardcoded
    transmission_rates::Vector{Float64}
    recovery_rates::Vector{Float64}
    imprinted_loss_rate::Float64
    matured_loss_rate::Float64
    frac_imprinted::Dict{Tuple{Vector{Int64},Int64},Float64}
    frac_matured::Dict{Tuple{Vector{Int64},Int64},Float64}
    immunity_coef::Dict{Tuple{Vector{Int64},Int64},Float64}
    evorisk_coef::Dict{Tuple{Vector{Int64},Int64},Float64}
    birth_rate_uninfected::Float64
    birth_rates_infected::Vector{Float64}
    death_rate_uninfected::Float64
    death_rates_infected::Vector{Float64}
end

function parameters(;
    n_strains::Int64,
    transmission_rates::Vector{Float64},
    recovery_rates::Vector{Float64},
    imprinted_loss_rate::Float64,
    matured_loss_rate::Float64,
    frac_imprinted::Matrix{Float64},
    frac_matured::Matrix{Float64},
    immunity_coef::Matrix{Float64},
    evorisk_coef::Array{Float64,3},
    birth_rate_uninfected::Float64,
    birth_rates_infected::Vector{Float64},
    death_rate_uninfected::Float64,
    death_rates_infected::Vector{Float64},)

    n_immunities = 3^n_strains # number of types of compartments, i.e. possible immune state combinations
    n_compartments = (n_strains + 1) * n_immunities + n_strains # n_strains possible infecting strains plus uninfected per each immune compartment type, plus one compartment per strain for evorisk motes

    immunities_ids = reduce(vcat, collect(multiset_permutations(i, n_strains)) for i in with_replacement_combinations([NAIVE, IMPRINTED, MATURED], n_strains))
    immunities_dict = Dict{Vector{Int64},Int64}(zip(immunities_ids, collect(1:n_immunities)))
    immunities_list = reduce(vcat, collect(multiset_permutations(i, n_strains)) for i in with_replacement_combinations([NAIVE_VEC, IMPRINTED_VEC, MATURED_VEC], n_strains))
    immunities = cat(dims=3, [hcat(compartment...) for compartment in immunities_list]...)
    # Dimensions of this matrix are: immune state (naive, imprinted, matured):strains:host compartment
    # This entire matrix appliers to both uninfected and infected with each strain, so there are (1 + n_strains) as many compartments as dimension 3 of this matrix

    immunity_coef_by_state = Dict{Tuple{Vector{Int64},Int64},Float64}()
    frac_imprinted_by_state = Dict{Tuple{Vector{Int64},Int64},Float64}()
    frac_matured_by_state = Dict{Tuple{Vector{Int64},Int64},Float64}()
    evorisk_by_state = Dict{Tuple{Vector{Int64},Int64},Float64}()
    for id in immunities_ids
        for s in 1:n_strains
            max_immunity_to_s = 0.0
            max_frac_matured = 0.0 # if no immunity, de novo immunity is imprinted
            corresp_frac_imprinted = 1.0 # if no immunity, de novo immunity is imprinted
            evorisk = evorisk_coef[s, s, 1] # if no immunity, evorisk is pulled from naive layer of array
            for i in 1:n_strains
                if id[i] > 0 && immunity_coef[i, s] > max_immunity_to_s
                    max_immunity_to_s = immunity_coef[i, s]
                    evorisk = evorisk_coef[i, s, id[i]+1]
                end
                if id[i] > 0 && frac_matured[i, s] > max_frac_matured
                    max_frac_matured = frac_matured[i, s]
                    corresp_frac_imprinted = frac_imprinted[i, s]
                end
            end
            immunity_coef_by_state[(id, s)] = 1 - max_immunity_to_s
            frac_matured_by_state[(id, s)] = max_frac_matured
            frac_imprinted_by_state[(id, s)] = corresp_frac_imprinted
            evorisk_by_state[(id, s)] = evorisk
        end
    end

    params = ModelParameters(
        n_strains=n_strains,
        n_immunities=n_immunities,
        n_compartments=n_compartments,
        immunities_ids=immunities_ids,
        immunities_dict=immunities_dict,
        immunities_list=immunities_list,
        immunities=immunities,
        transmission_rates=transmission_rates,
        recovery_rates=recovery_rates,
        imprinted_loss_rate=imprinted_loss_rate,
        matured_loss_rate=matured_loss_rate,
        frac_imprinted=frac_imprinted_by_state,
        frac_matured=frac_matured_by_state,
        immunity_coef=immunity_coef_by_state,
        evorisk_coef=evorisk_by_state,
        birth_rate_uninfected=birth_rate_uninfected,
        birth_rates_infected=birth_rates_infected,
        death_rate_uninfected=death_rate_uninfected,
        death_rates_infected=death_rates_infected,
    )

    return params
end
