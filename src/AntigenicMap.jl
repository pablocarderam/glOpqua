using Distances
using LinearAlgebra

@kwdef struct AntigenicMap
    positions::Matrix{Float64}
    max_cross_immunity::Float64
    K_cross_immunity::Float64
    h_cross_immunity::Float64
    max_imprinted::Float64
    K_imprinted::Float64
    h_imprinted::Float64
    max_matured::Float64
    K_matured::Float64
    h_matured::Float64

    distances::Matrix{Float64}
    # pairwise antigenic distances
    immunity_coef::Matrix{Float64}
    # This tracks strain cross immunity (self immunity is on diagonal).
    # Row number (first index) is pre-existing immunity to that number strain,
    # columns are number of infecting strain.
    frac_imprinted::Matrix{Float64}
    # This tracks likelihood of developing an imprinted response to a new strain given
    # pre-existing immunity to another strain.
    # Diagonal is zero (immunity is not re-imprinted on infection with same strain).
    # Row number (first index) is pre-existing immunity to that number strain,
    # columns are number of infecting strain.
    frac_matured::Matrix{Float64}
    # This tracks likelihood of developing an affinity matured response to a new strain given
    # pre-existing immunity to another strain.
    # Diagonal is zero (it is true immunity is affinity-matured on reexposure,
    # but negative effects on antigenic escape likelihood are avoided).
    # Row number (first index) is pre-existing immunity to that number strain,
    # columns are number of infecting strain.
    #
    # Notice frac_imprinted and frac_matured must sum to <=1 at each element.
    # 1-frac_imprinted-frac_matured is the fraction of cases with preexisting immunity
    # to the row number strain that do not acquire any immunity to the column number
    # strain after recovering from that infection (this diagonal is not taken into
    # account inside the ODE function).
end

function antigenicDistances(positions::Matrix{Float64}, distanceMetric=euclidean)
    distances = zeros(size(positions)[1], size(positions)[1])

    for s in 1:size(positions)[1]
        for s2 in 1:size(positions)[1]
            distances[s, s2] = distanceMetric(positions[s, :], positions[s2, :])
        end
    end

    return distances
end

function crossImmunity(max_cross_immunity::Float64, K_cross_immunity::Float64, h_cross_immunity::Float64, distances::Matrix{Float64})
    return max_cross_immunity * (1.0 .- hillFunction.(distances, K_cross_immunity, h_cross_immunity))
end

function imprinted(max_imprinted::Float64, K_imprinted::Float64, h_imprinted::Float64, distances::Matrix{Float64})
    return max_imprinted * hillFunction.(distances, K_imprinted, h_imprinted)
end

function matured(max_matured::Float64, K_matured::Float64, h_matured::Float64, distances::Matrix{Float64})
    return max_matured * (1.0 .- hillFunction.(distances, K_matured, h_matured)) .* (ones(size(distances)[1], size(distances)[2]) - Matrix{Float64}(I, size(distances)[1], size(distances)[2]))
end

function antigenicMap(;
    positions::Matrix{Float64},
    max_cross_immunity::Float64,
    K_cross_immunity::Float64,
    h_cross_immunity::Float64,
    max_imprinted::Float64,
    K_imprinted::Float64,
    h_imprinted::Float64,
    max_matured::Float64,
    K_matured::Float64,
    h_matured::Float64,)

    distances = antigenicDistances(positions)
    immunity_coef = crossImmunity(max_cross_immunity, K_cross_immunity, h_cross_immunity, distances)
    frac_imprinted = imprinted(max_imprinted, K_imprinted, h_imprinted, distances)
    frac_matured = matured(max_matured, K_matured, h_matured, distances)

    map = AntigenicMap(positions=positions,
        max_cross_immunity=max_cross_immunity,
        K_cross_immunity=K_cross_immunity,
        h_cross_immunity=h_cross_immunity,
        max_imprinted=max_imprinted,
        K_imprinted=K_imprinted,
        h_imprinted=h_imprinted,
        max_matured=max_matured,
        K_matured=K_matured,
        h_matured=h_matured,
        distances=distances,
        immunity_coef=immunity_coef,
        frac_imprinted=frac_imprinted,
        frac_matured=frac_matured,)

    return map
end
