function multistrain2Immunity(du, u, p::ModelParameters, t::Float64)
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
            infected[s] = p.transmission_rates[s] * p.immunity_coef[(p.immunities_ids[c],s)] * u[c] * sum(@views u[(s*p.n_immunities+1):((s+1)*p.n_immunities)])
            total_infected += infected[s]
            # It should be true that immunities_dict[immunities_ids[c]] == c, so instead of u[immunities_dict[immunities_ids[c]]+(s*n_immunities)] in these next three blocks, we do u[c+(s*n_immunities)]
            if p.immunities_ids[c][s] == IMPRINTED
                immunity = copy(p.immunities_ids[c])
                immunity[s] = 0 # We consider the immune status identical to the current one but with naive in the current strain
                imprinted += p.recovery_rates[s] * (p.frac_imprinted[(p.immunities_ids[c],s)] * u[p.immunities_dict[immunity]+(s*p.n_immunities)] + u[c+(s*p.n_immunities)])
                # imprinted += p.recovery_rates[s] * (0 * u[p.immunities_dict[immunity]+(s*p.n_immunities)] + u[c+(s*p.n_immunities)])
            elseif p.immunities_ids[c][s] == MATURED
                immunity = copy(p.immunities_ids[c])
                immunity[s] = 0 # We consider the immune status identical to the current one but with naive in the current strain
                matured += p.recovery_rates[s] * (p.frac_matured[(p.immunities_ids[c],s)] * u[p.immunities_dict[immunity]+(s*p.n_immunities)] + u[c+(s*p.n_immunities)])
                # matured += p.recovery_rates[s] * (1 * u[p.immunities_dict[immunity]+(s*p.n_immunities)] + u[c+(s*p.n_immunities)])
            elseif p.immunities_ids[c][s] == NAIVE
                naive += p.recovery_rates[s] * (1 - p.frac_imprinted[(p.immunities_ids[c],s)] - p.frac_matured[(p.immunities_ids[c],s)]) * u[c+(s*p.n_immunities)]
                #TODO: add loss of immunity (different rates for imprinted and matured), add population turnover
            end
        end

        # Uninfected compartment for this immunity type:
        du[c] = -total_infected + imprinted + matured + naive

        # Infected compartments for this immunity type:
        for s in 1:p.n_strains
            du[c+(s*p.n_immunities)] = infected[s] - p.recovery_rates[s] * u[c+(s*p.n_immunities)]
        end
    end
end
