using Combinatorics

const N_IMMUNE_STATES = 3

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
    n_ag_clusters::Int64
    n_strains::Vector{Int64}
    n_immunities::Int64
    n_compartments::Int64
    immunities_ids::Vector{Vector{Int64}}
    immunities_dict::Dict{Vector{Int64},Int64}
    immunities_list::Vector{Vector{Vector{Bool}}}
    # immunities::Array{Bool,N_IMMUNE_STATES}
    transmission_rates::Vector{Vector{Float64}}
    recovery_rates::Vector{Vector{Float64}}
    imprinted_loss_rate::Float64
    matured_loss_rate::Float64
    frac_imprinted::Dict{Tuple{Vector{Int64},Int64},Float64}
    frac_matured::Dict{Tuple{Vector{Int64},Int64},Float64}
    immunity_coef::Dict{Tuple{Vector{Int64},Int64},Float64}
    evorisk_coef::Dict{Tuple{Vector{Int64},Int64},Float64}
    birth_rate_uninfected::Float64
    birth_rates_infected::Vector{Vector{Float64}}
    death_rate_uninfected::Float64
    death_rates_infected::Vector{Vector{Float64}}
end

function parameters(;
    n_ag_clusters::Int64,
    n_strains::Vector{Int64},
    transmission_rates::Vector{Vector{Float64}},
    recovery_rates::Vector{Vector{Float64}},
    imprinted_loss_rate::Float64,
    matured_loss_rate::Float64,
    frac_imprinted::Matrix{Float64},
    frac_matured::Matrix{Float64},
    immunity_coef::Matrix{Float64},
    evorisk_coef::Array{Float64,3},
    birth_rate_uninfected::Float64,
    birth_rates_infected::Vector{Vector{Float64}},
    death_rate_uninfected::Float64,
    death_rates_infected::Vector{Vector{Float64}},)

    n_immunities = 3^n_ag_clusters # number of types of compartments, i.e. possible immune state combinations
    total_n_strains = sum(n_strains)
    n_compartments = (total_n_strains + 1) * n_immunities + total_n_strains # n_strains possible infecting strains plus uninfected per each immune compartment type, plus one compartment per strain for evorisk motes

    immunities_ids = reduce(vcat, collect(multiset_permutations(i, n_ag_clusters)) for i in with_replacement_combinations([NAIVE, IMPRINTED, MATURED], n_ag_clusters))
    immunities_dict = Dict{Vector{Int64},Int64}(zip(immunities_ids, collect(1:n_immunities)))
    immunities_list = reduce(vcat, collect(multiset_permutations(i, n_ag_clusters)) for i in with_replacement_combinations([NAIVE_VEC, IMPRINTED_VEC, MATURED_VEC], n_ag_clusters))
    println(length(immunities_list))
    # immunities = cat(dims=3, [hcat(compartment...) for compartment in immunities_list]...)
    println("Out")
    # Dimensions of this matrix are: immune state (naive, imprinted, matured):strains:host compartment
    # This entire matrix appliers to both uninfected and infected with each strain, so there are (1 + n_strains) as many compartments as dimension 3 of this matrix

    immunity_coef_by_state = Dict{Tuple{Vector{Int64},Int64},Float64}()
    frac_imprinted_by_state = Dict{Tuple{Vector{Int64},Int64},Float64}()
    frac_matured_by_state = Dict{Tuple{Vector{Int64},Int64},Float64}()
    evorisk_by_state = Dict{Tuple{Vector{Int64},Int64},Float64}()
    ag_cluster_s = 1
    cluster_counter_s = 0
    ag_cluster_i = 1
    cluster_counter_i = 0
    for id in immunities_ids
        for s in 1:n_ag_clusters
            max_immunity_to_s = 0.0
            max_frac_matured = 0.0 # if no immunity, de novo immunity is imprinted
            corresp_frac_imprinted = 1.0 # if no immunity, de novo immunity is imprinted
            for i in 1:n_ag_clusters
                if id[i] > 0 && immunity_coef[i, s] > max_immunity_to_s
                    max_immunity_to_s = immunity_coef[i, s]
                end
                if id[i] > 0 && frac_matured[i, s] > max_frac_matured
                    max_frac_matured = frac_matured[i, s]
                    corresp_frac_imprinted = frac_imprinted[i, s]
                end
            end
            immunity_coef_by_state[(id, s)] = 1 - max_immunity_to_s
            frac_matured_by_state[(id, s)] = max_frac_matured
            frac_imprinted_by_state[(id, s)] = corresp_frac_imprinted
        end

        ag_cluster_s = 1
        cluster_counter_s = 0
        for s in 1:total_n_strains
            cluster_counter_s += 1
            if cluster_counter_s > n_strains[ag_cluster_s]
                ag_cluster_s += 1
                cluster_counter_s = 0
            end
            max_immunity_to_s = 0.0
            evorisk = evorisk_coef[s, s, 1] # if no immunity to this strain's antigenic cluster, evorisk is pulled from naive layer of array
            ag_cluster_i = 1
            cluster_counter_i = 0
            for i in 1:total_n_strains
                cluster_counter_i += 1
                if cluster_counter_i > n_strains[ag_cluster_i]
                    ag_cluster_i += 1
                    cluster_counter_i = 0
                elseif id[ag_cluster_i] > 0 && immunity_coef[ag_cluster_i, ag_cluster_s] > max_immunity_to_s
                    max_immunity_to_s = immunity_coef[ag_cluster_i, ag_cluster_s]
                    evorisk = evorisk_coef[i, s, id[ag_cluster_s]+1]
                    # evorisk is determined by strongest immune reaction to this strain's antigenic cluster
                end
            end
            evorisk_by_state[(id, s)] = evorisk
        end
    end

    params = ModelParameters(
        n_ag_clusters=n_ag_clusters,
        n_strains=n_strains,
        n_immunities=n_immunities,
        n_compartments=n_compartments,
        immunities_ids=immunities_ids,
        immunities_dict=immunities_dict,
        immunities_list=immunities_list,
        # immunities=immunities,
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
