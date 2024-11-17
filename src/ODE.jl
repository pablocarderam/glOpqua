function multistrainImprintMature(du, u, p::ModelParameters, t::Float64)
    # Temp variables for the loop:
    infected::Vector{Float64} = zeros(Float64, p.n_immunities)
    immunity_naive::Vector{Int64} = zeros(Int64, p.n_strains)
    immunity_imprinted::Vector{Int64} = zeros(Int64, p.n_strains)
    immunity_matured::Vector{Int64} = zeros(Int64, p.n_strains)

    # Loop through all immune types of compartments
    for c in 1:p.n_immunities
        naive_uninfected = 0.0
        imprinted_uninfected = 0.0
        matured_uninfected = 0.0
        naive_infected = 0.0
        imprinted_infected = 0.0
        matured_infected = 0.0
        total_infected = 0.0
        # Calculate totals:
        for s in 1:p.n_strains
            infected[s] = p.transmission_rates[s] * p.immunity_coef[(p.immunities_ids[c], s)] * u[c] * sum(@views u[(s*p.n_immunities+1):((s+1)*p.n_immunities)])
            total_infected += infected[s]

            # It should be true that immunities_dict[immunities_ids[c]] == c, so instead of u[immunities_dict[immunities_ids[c]]+(s*n_immunities)] in these next three blocks, we do u[c+(s*n_immunities)]
            if p.immunities_ids[c][s] == IMPRINTED
                immunity_naive = copy(p.immunities_ids[c])
                immunity_naive[s] = NAIVE # We consider the immune status identical to the current one but with naive in the current strain
                imprinted_uninfected += -p.imprinted_loss_rate * u[c] + p.recovery_rates[s] * (p.frac_imprinted[(p.immunities_ids[c], s)] * u[p.immunities_dict[immunity_naive]+(s*p.n_immunities)] + u[c+(s*p.n_immunities)])

                imprinted_infected += -p.imprinted_loss_rate * u[c+(s*p.n_immunities)]
            elseif p.immunities_ids[c][s] == MATURED
                immunity_naive = copy(p.immunities_ids[c])
                immunity_naive[s] = NAIVE # We consider the immune status identical to the current one but with naive in the current strain
                matured_uninfected += -p.matured_loss_rate * u[c] + p.recovery_rates[s] * (p.frac_matured[(p.immunities_ids[c], s)] * u[p.immunities_dict[immunity_naive]+(s*p.n_immunities)] + u[c+(s*p.n_immunities)])

                matured_infected += -p.matured_loss_rate * u[c+(s*p.n_immunities)]
            elseif p.immunities_ids[c][s] == NAIVE
                immunity_imprinted = copy(p.immunities_ids[c])
                immunity_imprinted[s] = IMPRINTED # We consider the immune status identical to the current one but with imprinted in the current strain
                immunity_matured = copy(p.immunities_ids[c])
                immunity_matured[s] = MATURED # We consider the immune status identical to the current one but with naive in the current strain
                naive_uninfected += p.imprinted_loss_rate * u[p.immunities_dict[immunity_imprinted]] + p.imprinted_loss_rate * u[p.immunities_dict[immunity_matured]] + p.recovery_rates[s] * (1 - p.frac_imprinted[(p.immunities_ids[c], s)] - p.frac_matured[(p.immunities_ids[c], s)]) * u[c+(s*p.n_immunities)]

                naive_infected += p.imprinted_loss_rate * u[p.immunities_dict[immunity_imprinted]+(s*p.n_immunities)] + p.imprinted_loss_rate * u[p.immunities_dict[immunity_matured]+(s*p.n_immunities)]
            end
        end

        # Uninfected compartment for this immunity type:
        du[c] = -total_infected + imprinted_uninfected + matured_uninfected + naive_uninfected - p.death_rate_uninfected * u[c]

        # Add birth rate to fully naive, uninfected compartment:
        du[1] = du[1] + p.birth_rate_uninfected * u[c]

        # Infected compartments for this immunity type:
        for s in 1:p.n_strains
            du[c+(s*p.n_immunities)] = infected[s] - p.recovery_rates[s] * u[c+(s*p.n_immunities)] + imprinted_infected + matured_infected + naive_infected - p.death_rates_infected[s] * u[c+(s*p.n_immunities)]
            # Add birth rate to fully naive, uninfected compartment:
            du[1] = du[1] + p.birth_rates_infected[s] * u[c+(s*p.n_immunities)]
        end
    end


end
