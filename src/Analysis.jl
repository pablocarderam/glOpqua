using DataFrames

function dataByStrain(ode_sol::SciMLBase.ODESolution,params::ModelParameters)
    dat = DataFrame(Tables.table(ode_sol'))
    dat[!, :t] = ode_sol.t

    dat_agg = DataFrame()
    dat_agg[!, :t] = ode_sol.t

    dat_agg[!, ("Uninfected")] = sum(eachcol(dat[:, 1:params.n_immunities]))
    for i in 1:params.n_strains
        dat_agg[!, ("Strain_"*string(i))] = sum(eachcol(dat[:, (params.n_immunities*i+1):(params.n_immunities*(i+1))]))
    end

    return dat_agg
end
