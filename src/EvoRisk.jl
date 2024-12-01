using Distances
using LinearAlgebra

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

    evorisk_coef::Array{Float64,3}
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
    strains_per_cluster::Vector{Int64},)

    evorisk_coef = evorisk(evorisk_by_immunity, K_evorisk, h_evorisk, ag_distances, strains_per_cluster)

    map = EvoRiskMap(#positions=positions,
        K_evorisk=K_evorisk,
        h_evorisk=h_evorisk,
        evorisk_by_immunity=evorisk_by_immunity,
        ag_distances=ag_distances,
        evorisk_coef=evorisk_coef,)

    return map
end
