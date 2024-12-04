using Distances
using LinearAlgebra
using Random

@kwdef struct EvoRiskMap
    K_evorisk::Vector{Float64}
    h_evorisk::Vector{Float64}
    ag_distances::Matrix{Float64}
    # The following evorisks are related to each other and are functions of:
    # 1. the landscape surrounding a strain in given immune background (i.e. number and fitness of paths out)
    # 2. relative fitness of current strain in given immune background
    # 3. inoculum bottleneck effects
    # 4. stochastic extinction/survival probability before getting to deterministic epidemiological dynamic threshold
    evorisk_by_immunity::Matrix{Float64}
    strains_per_cluster::Vector{Int64}

    evorisk_coef::Array{Float64,3}
    #TODO: we need a separate evorisk_coef array for (non-cluster) strain evolution, or to specify the breakdown of cluster vs. non-cluster
end

function evorisk(evorisk_by_immunity::Matrix{Float64}, K_evorisk::Vector{Float64}, h_evorisk::Vector{Float64}, ag_distances::Matrix{Float64}, strains_per_cluster::Vector{Int64})
    evorisk_coef = zeros(length(strains_per_cluster), sum(strains_per_cluster), 3)
    strains_in_prev_clusters = 0
    for c in 1:length(strains_per_cluster)
        cluster_s = 1
        strain_count = 0
        for s in 1:sum(strains_per_cluster)
            strain_count += 1
            if strain_count > strains_per_cluster[cluster_s]
                cluster_s += 1
                strain_count = 1
            end
            evorisk_coef[c, s, NAIVE+1] = evorisk_by_immunity[NAIVE+1, s]
            evorisk_coef[c, s, IMPRINTED+1] = evorisk_by_immunity[IMPRINTED+1, s] * (1.0 - hillFunction(ag_distances[cluster_s, c], K_evorisk[IMPRINTED+1], h_evorisk[IMPRINTED+1]))
            evorisk_coef[c, s, MATURED+1] = evorisk_by_immunity[MATURED+1, s] * (1.0 - hillFunction(ag_distances[cluster_s, c], K_evorisk[MATURED+1], h_evorisk[MATURED+1]))
        end
        strains_in_prev_clusters += strains_per_cluster[c]
    end

    # # These notes below are not necessarily true and are just for record of previous thought processes:
    # # If there is 1 strain per cluster, this should be a diagonal matrix, one value per strain (assuming no evolution in naive hosts)
    # # In this antigenic, non-sequence based model, this is all zeros

    # # NOTE: In this antigenic, non-sequence based model, the diagonal of each of these is all the same and everything outside
    # # of the diagonal increases following a function that tracks the mass of the immune hypersphere formed
    # # by the mutant cloud of the infecting strain that lies inside of the volume of the hypersphere corresponding
    # # to the pre-existing cross-reactive immunity, up until the point in which the infecting strain itself falls
    # # outside the pre-existing immunity hypersphere, at which point we assume evorisk goes to 0.

    # #TODO: The function defining the modifying coefficients as a function of distance depends on
    # # simulations to be figured out in stochastic strain model and popgen+sim to be done in Opqua/jOpqua 2.0.
    # # For the purposes of this model, we assume antigenic clusters are separate enough from each other that anything outside the diagonal tends to zero.
    # # We can use a Hill function to define that.
    # evorisk_coef[:, :, IMPRINTED+1] = imprinted_matured_evorisk[IMPRINTED] * (1.0 .- hillFunction.(ag_distances, K_evorisk[IMPRINTED+1], h_evorisk[IMPRINTED+1]))
    # # This tracks evolutionary risk motes for imprinted immunity.
    # # Row number (first index) is pre-existing immunity to that number strain,
    # # columns are number of infecting strain.
    # evorisk_coef[:, :, MATURED+1] = imprinted_matured_evorisk[MATURED] * (1.0 .- hillFunction.(ag_distances, K_evorisk[MATURED+1], h_evorisk[MATURED+1]))
    # # This tracks evolutionary risk motes for affinity matured immunity.
    # # Row number (first index) is pre-existing immunity to that number strain,
    # # columns are number of infecting strain.

    return evorisk_coef
end

function evoRiskMap(;
    ag_distances::Matrix{Float64},
    K_evorisk::Vector{Float64},
    h_evorisk::Vector{Float64},
    evorisk_by_immunity::Matrix{Float64},
    strains_per_cluster::Vector{Int64},
)

    evorisk_coef = evorisk(evorisk_by_immunity, K_evorisk, h_evorisk, ag_distances, strains_per_cluster)

    map = EvoRiskMap(#positions=positions,
        K_evorisk=K_evorisk,
        h_evorisk=h_evorisk,
        evorisk_by_immunity=evorisk_by_immunity,
        ag_distances=ag_distances,
        evorisk_coef=evorisk_coef,
        strains_per_cluster=strains_per_cluster,
    )

    return map
end

function copyEvoRiskMap(p::EvoRiskMap)
    return evoRiskMap(
        ag_distances=p.ag_distances,
        K_evorisk=p.K_evorisk,
        h_evorisk=p.h_evorisk,
        evorisk_by_immunity=p.evorisk_by_immunity,
        strains_per_cluster=p.strains_per_cluster,
    )
end

function sampleInfection(sol, p::ModelParameters, ag::AntigenicMap)
    u = sol.u[end]
    rates = zeros(size(sol.u[end])[1] * 2)#p.n_immunities * sum(p.n_strains) * 2)
    ids = [(zeros(p.n_ag_clusters), zeros(size(ag.positions)[2]), 0, 0, false) for _ in rates]
    # Loop through all immune types of compartments
    for c in 1:p.n_immunities
        strain_num = 0
        for a in 1:p.n_ag_clusters
            for s in 1:p.n_strains[a]
                strain_num += 1
                # Evolutionary risk motes:
                # New cluster:
                rates[c+(strain_num*p.n_immunities)] = p.transmission_rates[a][s] * p.evorisk_coef[(p.immunities_ids[c], strain_num)] * u[c] * sum(@views u[(strain_num*p.n_immunities+1):((strain_num+1)*p.n_immunities)])
                ids[c+(strain_num*p.n_immunities)] = (p.immunities_ids[c], ag.positions[a,:], a, s, false)
                # New strain in old cluster:
                # println((size(sol.u[end]), strain_num, c + (strain_num * p.n_immunities), c + (strain_num * p.n_immunities) + (p.n_immunities * sum(p.n_strains))))
                rates[c+(strain_num*p.n_immunities)+(p.n_immunities*sum(p.n_strains))] = 0.0*p.transmission_rates[a][s] * p.evorisk_coef[(p.immunities_ids[c], strain_num)] * u[c] * sum(@views u[(strain_num*p.n_immunities+1):((strain_num+1)*p.n_immunities)])
                #TODO: this is currently set to zero but needs to accommodate new strain evo
                ids[c+(strain_num*p.n_immunities)+(p.n_immunities*sum(p.n_strains))] = (p.immunities_ids[c], ag.positions[a,:], a, s, true)
            end
        end
    end

    id = sample(ids, Weights(rates), 1)[1]
    infecting_host_immunities_id = id[1]

    parent_cluster = id[2]
    parent_cluster_num = id[3]
    parent_strain = id[4]
    same_clusters = id[5]

    infecting_host_immunities = [zeros(size(ag.positions)[2]) for _ in 1:sum(id[1] .> NAIVE)]
    counter = 0
    for i in 1:p.n_ag_clusters
        if id[1][i] > 0
            counter += 1
            infecting_host_immunities[counter] = ag.positions[i, :]
        end
    end

    return parent_strain, parent_cluster, parent_cluster_num, infecting_host_immunities, infecting_host_immunities_id, same_clusters
end

function mutantFromParent(parent_num::Int64, parent_cluster::Vector{Float64}, parent_cluster_num::Int64, immunities::Vector{Vector{Float64}}, same_clusters::Bool; distanceFunc::Function=(p) -> randexp(), distanceMetric=euclidean, max_attempts::Int64=100)
    if same_clusters
        return mutantStrainFromParent(), parent_cluster_num
    else
        return mutantAgFromParent(parent_cluster, immunities; distanceFunc=distanceFunc, distanceMetric=distanceMetric, max_attempts=max_attempts),-1
    end
end

function mutantStrainFromParent()
    return 0
end

function mutantAgFromParent(parent::Vector{Float64}, other_strains::Vector{Vector{Float64}}; distanceFunc::Function=(p) -> randexp(), distanceMetric=euclidean, max_attempts::Int64=100)
    parent_distances = [distanceMetric(parent, s) for s in other_strains]
    mutant = mutantAgFromParent(parent, distanceFunc=distanceFunc)
    attempts = 0
    println(parent_distances)
    println(([distanceMetric(mutant, s) for s in other_strains], [distanceMetric(mutant, s) for s in other_strains] .< parent_distances))
    while any([distanceMetric(mutant, s) for s in other_strains] .< parent_distances) && attempts < max_attempts
        attempts += 1
        println((attempts, [distanceMetric(mutant, s) for s in other_strains], [distanceMetric(mutant, s) for s in other_strains] .< parent_distances))
        mutant = mutantAgFromParent(parent, distanceFunc=distanceFunc)
    end

    return mutant
end

function mutantAgFromParent(parent::Vector{Float64}; distanceFunc::Function=(p) -> randexp())
    dimensions = length(parent)
    angles = pi * rand(dimensions - 1)
    sines = sin.(angles)
    cosines = cos.(angles)
    norm_coordinates = [cosines[1], [prod(sines[1:i-1]) * cosines[i] for i in 2:dimensions-1]..., prod(sines)]
    mutation_vec = randexp() .* norm_coordinates

    return parent .+ mutation_vec
end
