using Distances
using LinearAlgebra

@kwdef struct EvoRiskMap
    positions::Matrix{Float64}
    K_evorisk::Vector{Float64}
    h_evorisk::Vector{Float64}

    # Both of the following evorisks are related to each other and are functions of:
    # 1. the landscape surrounding a strain in given immune background (i.e. number and fitness of paths out)
    # 2. relative fitness of current strain in given immune background
    # 3. inoculum bottleneck effects
    # 4. stochastic extinction/survival probability before getting to deterministic epidemiological dynamic threshold
    naive_evorisk::Vector{Float64}
    imprinted_matured_evorisk::Matrix{Float64}

    distances::Matrix{Float64}
    evorisk_coef::Array{Float64,3}
end

function evoDistances(positions::Matrix{Float64}, distanceMetric=euclidean)
    distances = zeros(size(positions)[1], size(positions)[1])

    for s in 1:size(positions)[1]
        for s2 in 1:size(positions)[1]
            distances[s, s2] = distanceMetric(positions[s, :], positions[s2, :])
        end
    end

    return distances
end

function evorisk(naive_evorisk::Vector{Float64}, imprinted_matured_evorisk::Matrix{Float64}, K_evorisk::Vector{Float64}, h_evorisk::Vector{Float64}, distances::Matrix{Float64})
    evorisk_coef = zeros(size(distances)[1], size(distances)[1], 3)
    evorisk_coef[:, :, NAIVE+1] = diagm(naive_evorisk)
    # Should be a diagonal matrix, one value per strain
    # In this antigenic, non-sequence based model, this is all zeros

    # NOTE: In this antigenic, non-sequence based model, the diagonal of each of these is all the same and everything outside
    # of the diagonal increases following a function that tracks the mass of the immune hypersphere formed
    # by the mutant cloud of the infecting strain that lies inside of the volume of the hypersphere corresponding
    # to the pre-existing cross-reactive immunity, up until the point in which the infecting strain itself falls
    # outside the pre-existing immunity hypersphere, at which point we assume evorisk goes to 0.

    #TODO: The function defining the modifying coefficients as a function of distance depends on
    # simulations to be figured out in stochastic strain model and popgen+sim to be done in Opqua/jOpqua 2.0.
    # For the purposes of this model, we assume antigenic clusters are separate enough from each other that anything outside the diagonal tends to zero.
    # We can use a Hill function to define that.
    evorisk_coef[:, :, IMPRINTED+1] = imprinted_matured_evorisk[IMPRINTED] * (1.0 .- hillFunction.(distances, K_evorisk[2], h_evorisk[2]))
    # This tracks evolutionary risk motes for imprinted immunity.
    # Row number (first index) is pre-existing immunity to that number strain,
    # columns are number of infecting strain.
    evorisk_coef[:, :, MATURED+1] = imprinted_matured_evorisk[MATURED] * (1.0 .- hillFunction.(distances, K_evorisk[3], h_evorisk[3]))
    # This tracks evolutionary risk motes for affinity matured immunity.
    # Row number (first index) is pre-existing immunity to that number strain,
    # columns are number of infecting strain.

    return evorisk_coef
end

function evoRiskMap(;
    positions::Matrix{Float64},
    K_evorisk::Vector{Float64},
    h_evorisk::Vector{Float64},
    naive_evorisk::Vector{Float64},
    imprinted_matured_evorisk::Matrix{Float64},)

    distances = evoDistances(positions)
    evorisk_coef = evorisk(naive_evorisk, imprinted_matured_evorisk, K_evorisk, h_evorisk, distances)

    println(positions)
    println(K_evorisk)
    println(h_evorisk)
    println(naive_evorisk)
    println(imprinted_matured_evorisk)
    println(distances)
    println(evorisk_coef)

    map = EvoRiskMap(positions=positions,
        K_evorisk=K_evorisk,
        h_evorisk=h_evorisk,
        naive_evorisk=naive_evorisk,
        imprinted_matured_evorisk=imprinted_matured_evorisk,
        distances=distances,
        evorisk_coef=evorisk_coef,)

    return map
end
