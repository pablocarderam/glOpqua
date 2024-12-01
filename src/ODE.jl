function multistrainImprintMature(du, u, p::ModelParameters, t::Float64)
    # Temp variables for the loop:
    infected::Vector{Float64} = zeros(Float64, sum(p.n_strains))
    immunity_naive::Vector{Int64} = zeros(Int64, p.n_ag_clusters)
    immunity_imprinted::Vector{Int64} = zeros(Int64, p.n_ag_clusters)
    immunity_matured::Vector{Int64} = zeros(Int64, p.n_ag_clusters)

    # Reset these trackers
    for s in 1:sum(p.n_strains)
        du[(sum(p.n_strains)+1)*p.n_immunities+s] = 0.0
    end

    # Loop through all immune types of compartments
    for c in 1:p.n_immunities
        naive_uninfected = 0.0
        imprinted_uninfected = 0.0
        matured_uninfected = 0.0
        naive_infected_gained = 0.0
        imprinted_infected_lost = 0.0
        matured_infected_lost = 0.0
        total_infected = 0.0
        imprinted_lost = 0.0
        matured_lost = 0.0
        naive_gained = 0.0

        # Calculate totals:
        strain_num = 0
        for a in 1:p.n_ag_clusters
            for s in 1:p.n_strains[a]
                strain_num += 1
                infected[strain_num] = p.transmission_rates[a][s] * p.immunity_coef[(p.immunities_ids[c], a)] * u[c] * sum(@views u[(strain_num*p.n_immunities+1):((strain_num+1)*p.n_immunities)])
                total_infected += infected[strain_num]

                # It should be true that immunities_dict[immunities_ids[c]] == c, so instead of u[immunities_dict[immunities_ids[c]]+(s*n_immunities)] in these next three blocks, we do u[c+(s*n_immunities)]
                if p.immunities_ids[c][a] == IMPRINTED
                    immunity_naive = copy(p.immunities_ids[c])
                    immunity_naive[a] = NAIVE # We consider the immune status identical to the current one but with naive in the current strain
                    imprinted_uninfected += p.recovery_rates[a][s] * (p.frac_imprinted[(p.immunities_ids[c], a)] * u[p.immunities_dict[immunity_naive]+(strain_num*p.n_immunities)] + u[c+(strain_num*p.n_immunities)])

                    imprinted_lost += p.imprinted_loss_rate * (1.0 / p.n_strains[a]) * u[c]
                    imprinted_infected_lost += p.imprinted_loss_rate * u[c+(strain_num*p.n_immunities)]
                elseif p.immunities_ids[c][a] == MATURED
                    immunity_naive = copy(p.immunities_ids[c])
                    immunity_naive[a] = NAIVE # We consider the immune status identical to the current one but with naive in the current strain
                    matured_uninfected += p.recovery_rates[a][s] * (p.frac_matured[(p.immunities_ids[c], a)] * u[p.immunities_dict[immunity_naive]+(strain_num*p.n_immunities)] + u[c+(strain_num*p.n_immunities)])

                    matured_lost += p.matured_loss_rate * (1.0 / p.n_strains[a]) * u[c]
                    matured_infected_lost += p.matured_loss_rate * u[c+(strain_num*p.n_immunities)]
                elseif p.immunities_ids[c][a] == NAIVE
                    immunity_imprinted = copy(p.immunities_ids[c])
                    immunity_imprinted[a] = IMPRINTED # We consider the immune status identical to the current one but with imprinted in the current strain
                    immunity_matured = copy(p.immunities_ids[c])
                    immunity_matured[a] = MATURED # We consider the immune status identical to the current one but with naive in the current strain
                    naive_uninfected += p.recovery_rates[a][s] * (1 - p.frac_imprinted[(p.immunities_ids[c], a)] - p.frac_matured[(p.immunities_ids[c], a)]) * u[c+(strain_num*p.n_immunities)]

                    naive_gained += (1.0 / p.n_strains[a]) * (p.imprinted_loss_rate * u[p.immunities_dict[immunity_imprinted]] + p.matured_loss_rate * u[p.immunities_dict[immunity_matured]])
                    naive_infected_gained += p.imprinted_loss_rate * u[p.immunities_dict[immunity_imprinted]+(strain_num*p.n_immunities)] + p.matured_loss_rate * u[p.immunities_dict[immunity_matured]+(strain_num*p.n_immunities)]
                end
            end
        end

        # Uninfected compartment for this immunity type:
        du[c] = -total_infected + imprinted_uninfected + matured_uninfected + naive_uninfected - imprinted_lost - matured_lost + naive_gained - p.death_rate_uninfected * u[c]

        # Add birth rate to fully naive, uninfected compartment:
        du[1] = du[1] + p.birth_rate_uninfected * u[c]

        strain_num = 0
        for a in 1:p.n_ag_clusters
            for s in 1:p.n_strains[a]
                strain_num += 1
                # Infected compartments for this immunity type:
                du[c+(strain_num*p.n_immunities)] = infected[strain_num] - p.recovery_rates[a][s] * u[c+(strain_num*p.n_immunities)] - imprinted_infected_lost - matured_infected_lost + naive_infected_gained - p.death_rates_infected[a][s] * u[c+(strain_num*p.n_immunities)]
                # Add birth rate to fully naive, uninfected compartment:
                du[1] = du[1] + p.birth_rates_infected[a][s] * u[c+(strain_num*p.n_immunities)]

                # Evolutionary risk motes:
                du[(sum(p.n_strains)+1)*p.n_immunities+strain_num] += p.transmission_rates[a][s] * p.evorisk_coef[(p.immunities_ids[c], strain_num)] * u[c] * sum(@views u[(strain_num*p.n_immunities+1):((strain_num+1)*p.n_immunities)])
            end
        end
    end
end
