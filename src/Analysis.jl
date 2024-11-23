using DataFrames

function dataByStrain(ode_sol::SciMLBase.ODESolution, parameters::ModelParameters)
    dat = DataFrame(Tables.table(ode_sol'))
    dat[!, :t] = ode_sol.t

    dat_agg = DataFrame()
    dat_agg[!, :t] = ode_sol.t

    dat_agg[!, ("Uninfected")] = sum(eachcol(dat[:, 1:parameters.n_immunities]))
    for i in 1:parameters.n_strains
        dat_agg[!, ("Strain_"*string(i))] = sum(eachcol(dat[:, (parameters.n_immunities*i+1):(parameters.n_immunities*(i+1))]))
    end

    return dat_agg
end

function evoriskData(ode_sol::SciMLBase.ODESolution, parameters::ModelParameters)
    dat = DataFrame(Tables.table(ode_sol'))
    dat[!, :t] = ode_sol.t

    dat_evorisk = DataFrame()
    dat_evorisk[!, :t] = ode_sol.t

    for i in 1:parameters.n_strains
        dat_evorisk[!, ("Strain_"*string(i))] = dat[:, (parameters.n_strains+1)*parameters.n_immunities+i]
    end

    return dat_evorisk
end
